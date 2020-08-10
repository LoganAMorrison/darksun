//
// Created by logan on 8/9/20.
//

#ifndef DARKSUN_STIFF_RADAU_SOLVERS_HPP
#define DARKSUN_STIFF_RADAU_SOLVERS_HPP

#include "darksun/stiff/common.hpp"
#include "darksun/stiff/decomp.hpp"
#include "darksun/stiff/decsol.hpp"
#include "darksun/stiff/vector_matrix.hpp"
#include <cmath>

namespace darksun::stiff {

template <int N>
void slvrar(const Matrix<double, N, N> &fjac, const double fac1,
            const Matrix<double, N, N> &e1, Vector<double, N> &z1,
            const Vector<double, N> &f1, const Vector<int, N> &ip1) {
  for (int i = 0; i < N; ++i) {
    z1(i) -= f1(i) * fac1;
  }
  sol(e1, z1, ip1);
}
template <int N>
void slvrai(const Matrix<double, N, N> &fjac, const double alphn,
            const double betan, const Matrix<double, N, N> &e2r,
            const Matrix<double, N, N> &e2i, Vector<double, N> &z2,
            Vector<double, N> &z3, const Vector<double, N> &f2,
            const Vector<double, N> &f3, const Vector<int, N> &ip2) {
  for (int i = 0; i < N; ++i) {
    const double s2 = -f2(i);
    const double s3 = -f3(i);
    z2(i) = z2(i) + s2 * alphn - s3 * betan;
    z3(i) = z3(i) + s3 * alphn + s2 * betan;
  }
  solc(e2r, e2i, z2, z3, ip2);
}

template <int N>
void slvrad(const Matrix<double, N, N> &fjac, const double fac1,
            const double alphn, const double betan,
            const Matrix<double, N, N> &e1, const Matrix<double, N, N> &e2r,
            const Matrix<double, N, N> &e2i, Vector<double, N> &z1,
            Vector<double, N> &z2, Vector<double, N> &z3,
            const Vector<double, N> &f1, const Vector<double, N> &f2,
            const Vector<double, N> &f3, Vector<int, N> &ip1,
            Vector<int, N> &ip2) {
  for (int i = 0; i < N; ++i) {
    const double s2 = -f2(i);
    const double s3 = -f3(i);
    z1(i) -= f1(i) * fac1;
    z2(i) = z2(i) + s2 * alphn - s3 * betan;
    z3(i) = z3(i) + s3 * alphn + s2 * betan;
  }
  sol(e1, z1, ip1);
  solc(e2r, e2i, z2, z3, ip2);
}

} // namespace darksun::stiff

#endif // DARKSUN_STIFF_RADAU_SOLVERS_HPP
