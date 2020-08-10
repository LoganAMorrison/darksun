#ifndef DARKSUN_STIFF_DECSOL_HPP
#define DARKSUN_STIFF_DECSOL_HPP

#include "darksun/stiff/vector_matrix.hpp"

namespace darksun::stiff {

template <typename T, int SizeN = Eigen::Dynamic, int SizeM = Eigen::Dynamic>
int dec(Matrix<T, SizeN, SizeM> &a, Vector<int, SizeN> &ip) {
  /* System generated locals */
  int a_dim1, a_offset, i1, i2, i3;
  T d1, d2;

  /* Local variables */
  int i, j, k, m;
  T t;
  int nm1, kp1, ier;

  /* VERSION REAL DOUBLE PRECISION */
  /* ----------------------------------------------------------------------- */
  /*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION. */
  /*  INPUT.. */
  /*     N = ORDER OF MATRIX. */
  /*     NDIM = DECLARED DIMENSION OF ARRAY  A . */
  /*     A = MATRIX TO BE TRIANGULARIZED. */
  /*  OUTPUT.. */
  /*     A(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U . */
  /*     A(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L. */
  /*     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW. */
  /*     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O . */
  /*     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE */
  /*           SINGULAR AT STAGE K. */
  /*  USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM. */
  /*  DETERM(A) = IP(N)*A(1,1)*A(2,2)*...*A(N,N). */
  /*  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO. */

  /*  REFERENCE.. */
  /*     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, */
  /*     C.A.C.M. 15 (1972), P. 274. */
  /* ----------------------------------------------------------------------- */

  const int n = a.rows();
  const int ndim = a.cols();
  /* Function Body */
  ier = 0;
  ip(n - 1) = 1;
  if (n == 1) {
    goto L70;
  }
  nm1 = n - 1;
  i1 = nm1;
  for (k = 1; k <= i1; ++k) {
    kp1 = k + 1;
    m = k;
    i2 = n;
    for (i = kp1; i <= i2; ++i) {
      if ((d1 = a(i - 1, k - 1), abs(d1)) > (d2 = a(m - 1, k - 1), abs(d2))) {
        m = i;
      }
      /* L10: */
    }
    ip(k - 1) = m;
    t = a(m - 1, k - 1);
    if (m == k) {
      goto L20;
    }
    ip(n - 1) = -ip(n - 1);
    a(m - 1, k - 1) = a(k - 1, k - 1);
    a(k - 1, k - 1) = t;
  L20:
    if (t == 0.) {
      goto L80;
    }
    t = 1. / t;
    i2 = n;
    for (i = kp1; i <= i2; ++i) {
      /* L30: */
      a(i - 1, k - 1) = -a(i - 1, k - 1) * t;
    }
    i2 = n;
    for (j = kp1; j <= i2; ++j) {
      t = a(m - 1, j - 1);
      a(m - 1, j - 1) = a(k - 1, j - 1);
      a(k - 1, j - 1) = t;
      if (t == 0.) {
        goto L45;
      }
      i3 = n;
      for (i = kp1; i <= i3; ++i) {
        /* L40: */
        a(i - 1, j - 1) += a(i - 1, k - 1) * t;
      }
    L45:
        /* L50: */
        ;
    }
    /* L60: */
  }
L70:
  k = n;
  if (a(n - 1, n - 1) == 0.) {
    goto L80;
  }
  return ier;
L80:
  ier = k;
  ip(n - 1) = 0;
  return ier;
  /* ----------------------- END OF SUBROUTINE DEC ------------------------- */
}

template <typename T, int SizeN = Eigen::Dynamic, int SizeM = Eigen::Dynamic>
void sol(const Matrix<T, SizeN, SizeM> &a, Vector<T, SizeN> &b,
         const Vector<int, SizeN> &ip) {
  /* System generated locals */
  int a_dim1, a_offset, i1, i2;

  /* Local variables */
  int i, k, m;
  T t;
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
  const int n = a.rows();
  const int ndim = a.cols();

  /* Function Body */
  if (n == 1) {
    goto L50;
  }
  nm1 = n - 1;
  i1 = nm1;
  for (k = 1; k <= i1; ++k) {
    kp1 = k + 1;
    m = ip(k - 1);
    t = b(m - 1);
    b(m - 1) = b(k - 1);
    b(k - 1) = t;
    i2 = n;
    for (i = kp1; i <= i2; ++i) {
      /* L10: */
      b(i - 1) += a(i - 1, k - 1) * t;
    }
    /* L20: */
  }
  i1 = nm1;
  for (kb = 1; kb <= i1; ++kb) {
    km1 = n - kb;
    k = km1 + 1;
    b(k - 1) /= a(k - 1, k - 1);
    t = -b(k - 1);
    i2 = km1;
    for (i = 1; i <= i2; ++i) {
      /* L30: */
      b(i - 1) += a(i - 1, k - 1) * t;
    }
    /* L40: */
  }
L50:
  b(1 - 1) /= a(1 - 1, 1 - 1);
}

template <typename T, int SizeN = Eigen::Dynamic, int SizeM = Eigen::Dynamic>
int decc(Matrix<T, SizeN, SizeM> &ar, Matrix<T, SizeN, SizeM> &ai,
         Vector<int, SizeN> &ip) {
  /* System generated locals */
  int ar_dim1, ar_offset, ai_dim1, ai_offset, i1, i2, i3;
  T d1, d2, d__3, d__4;

  /* Local variables */
  int i, j, k, m;
  T ti, tr;
  int nm1, kp1, ier;
  T den, prodi, prodr;

  /* VERSION COMPLEX DOUBLE PRECISION */
  /* ----------------------------------------------------------------------- */
  /*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION */
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
  const int n = ar.rows();
  const int ndim = ar.cols();

  /* Function Body */
  ier = 0;
  ip(n - 1) = 1;
  if (n == 1) {
    goto L70;
  }
  nm1 = n - 1;
  i1 = nm1;
  for (k = 1; k <= i1; ++k) {
    kp1 = k + 1;
    m = k;
    i2 = n;
    for (i = kp1; i <= i2; ++i) {
      if ((d1 = ar(i - 1, k - 1), abs(d1)) + (d2 = ai(i - 1, k - 1), abs(d2)) >
          (d__3 = ar(m - 1, k - 1), abs(d__3)) +
              (d__4 = ai(m - 1, k - 1), abs(d__4))) {
        m = i;
      }
      /* L10: */
    }
    ip(k - 1) = m;
    tr = ar(m - 1, k - 1);
    ti = ai(m - 1, k - 1);
    if (m == k) {
      goto L20;
    }
    ip(n - 1) = -ip(n - 1);
    ar(m - 1, k - 1) = ar(k - 1, k - 1);
    ai(m - 1, k - 1) = ai(k - 1, k - 1);
    ar(k - 1, k - 1) = tr;
    ai(k - 1, k - 1) = ti;
  L20:
    if (abs(tr) + abs(ti) == 0.) {
      goto L80;
    }
    den = tr * tr + ti * ti;
    tr /= den;
    ti = -ti / den;
    i2 = n;
    for (i = kp1; i <= i2; ++i) {
      prodr = ar(i - 1, k - 1) * tr - ai(i - 1, k - 1) * ti;
      prodi = ai(i - 1, k - 1) * tr + ar(i - 1, k - 1) * ti;
      ar(i - 1, k - 1) = -prodr;
      ai(i - 1, k - 1) = -prodi;
      /* L30: */
    }
    i2 = n;
    for (j = kp1; j <= i2; ++j) {
      tr = ar(m - 1, j - 1);
      ti = ai(m - 1, j - 1);
      ar(m - 1, j - 1) = ar(k - 1, j - 1);
      ai(m - 1, j - 1) = ai(k - 1, j - 1);
      ar(k - 1, j - 1) = tr;
      ai(k - 1, j - 1) = ti;
      if (abs(tr) + abs(ti) == 0.) {
        goto L48;
      }
      if (ti == 0.) {
        i3 = n;
        for (i = kp1; i <= i3; ++i) {
          prodr = ar(i - 1, k - 1) * tr;
          prodi = ai(i - 1, k - 1) * tr;
          ar(i - 1, j - 1) += prodr;
          ai(i - 1, j - 1) += prodi;
          /* L40: */
        }
        goto L48;
      }
      if (tr == 0.) {
        i3 = n;
        for (i = kp1; i <= i3; ++i) {
          prodr = -ai(i - 1, k - 1) * ti;
          prodi = ar(i - 1, k - 1) * ti;
          ar(i - 1, j - 1) += prodr;
          ai(i - 1, j - 1) += prodi;
          /* L45: */
        }
        goto L48;
      }
      i3 = n;
      for (i = kp1; i <= i3; ++i) {
        prodr = ar(i - 1, k - 1) * tr - ai(i - 1, k - 1) * ti;
        prodi = ai(i - 1, k - 1) * tr + ar(i - 1, k - 1) * ti;
        ar(i - 1, j - 1) += prodr;
        ai(i - 1, j - 1) += prodi;
        /* L47: */
      }
    L48:
        /* L50: */
        ;
    }
    /* L60: */
  }
L70:
  k = n;
  if ((d1 = ar(n - 1, n - 1), abs(d1)) + (d2 = ai(n - 1, n - 1), abs(d2)) ==
      0.) {
    goto L80;
  }
  return ier;
L80:
  ier = k;
  ip(n - 1) = 0;
  return ier;
}

template <typename T, int SizeN = Eigen::Dynamic, int SizeM = Eigen::Dynamic>
void solc(const Matrix<T, SizeN, SizeM> &ar, const Matrix<T, SizeN, SizeM> &ai,
          Vector<T, SizeN> &br, Vector<T, SizeN> &bi,
          const Vector<int, SizeN> &ip) {
  /* System generated locals */
  int ar_dim1, ar_offset, ai_dim1, ai_offset, i1, i2;

  /* Local variables */
  int i, k, m, kb;
  T ti, tr;
  int km1, nm1, kp1;
  T den, prodi, prodr;

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
  const int n = ar.rows();

  /* Function Body */
  if (n == 1) {
    goto L50;
  }
  nm1 = n - 1;
  i1 = nm1;
  for (k = 1; k <= i1; ++k) {
    kp1 = k + 1;
    m = ip(k - 1);
    tr = br(m - 1);
    ti = bi(m - 1);
    br(m - 1) = br(k - 1);
    bi(m - 1) = bi(k - 1);
    br(k - 1) = tr;
    bi(k - 1) = ti;
    i2 = n;
    for (i = kp1; i <= i2; ++i) {
      prodr = ar(i - 1, k - 1) * tr - ai(i - 1, k - 1) * ti;
      prodi = ai(i - 1, k - 1) * tr + ar(i - 1, k - 1) * ti;
      br(i - 1) += prodr;
      bi(i - 1) += prodi;
      /* L10: */
    }
    /* L20: */
  }
  i1 = nm1;
  for (kb = 1; kb <= i1; ++kb) {
    km1 = n - kb;
    k = km1 + 1;
    den = ar(k - 1, k - 1) * ar(k - 1, k - 1) +
          ai(k - 1, k - 1) * ai(k - 1, k - 1);
    prodr = br(k - 1) * ar(k - 1, k - 1) + bi(k - 1) * ai(k - 1, k - 1);
    prodi = bi(k - 1) * ar(k - 1, k - 1) - br(k - 1) * ai(k - 1, k - 1);
    br(k - 1) = prodr / den;
    bi(k - 1) = prodi / den;
    tr = -br(k - 1);
    ti = -bi(k - 1);
    i2 = km1;
    for (i = 1; i <= i2; ++i) {
      prodr = ar(i - 1, k - 1) * tr - ai(i - 1, k - 1) * ti;
      prodi = ai(i - 1, k - 1) * tr + ar(i - 1, k - 1) * ti;
      br(i - 1) += prodr;
      bi(i - 1) += prodi;
      /* L30: */
    }
    /* L40: */
  }
L50:
  den = ar(0, 0) * ar(0, 0) + ai(0, 0) * ai(0, 0);
  prodr = br(0) * ar(0, 0) + bi(0) * ai(0, 0);
  prodi = bi(0) * ar(0, 0) - br(0) * ai(0, 0);
  br(0) = prodr / den;
  bi(0) = prodi / den;
}

} // namespace darksun::stiff

#endif // DARKSUN_STIFF_DECSOL_HPP
