#ifndef DARKSUN_STIFF_DECOMP_HPP
#define DARKSUN_STIFF_DECOMP_HPP

#include "darksun/stiff/common.hpp"
#include "darksun/stiff/decsol.hpp"
#include "darksun/stiff/vector_matrix.hpp"

namespace darksun::stiff {
template <int N>
int decomr(const Matrix<double, N, N> &fjac, const double fac1,
           Matrix<double, N, N> &e1, Vector<int, N> &ip1) {

  for (int j = 1; j <= N; ++j) {
    for (int i = 1; i <= N; ++i) {
      e1(i - 1, j - 1) = -fjac(i - 1, j - 1);
    }
    e1(j - 1, j - 1) += fac1;
  }
  return dec(e1, ip1);
}
template <int N>
int decomc(const Matrix<double, N, N> &fjac, const double alphn,
           const double betan, Matrix<double, N, N> &e2r,
           Matrix<double, N, N> &e2i, Vector<int, N> &ip2) {
  for (int j = 0; j < N; ++j) {
    for (int i = 0; i < N; ++i) {
      e2r(i, j) = -fjac(i, j);
      e2i(i, j) = 0.0;
    }
    e2r(j, j) += alphn;
    e2i(j, j) = betan;
  }
  return decc(e2r, e2i, ip2);
}
} // namespace darksun::stiff

#endif // DARKSUN_STIFF_DECOMP_HPP
