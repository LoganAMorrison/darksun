#ifndef DARKSUN_STIFF_RADAU_ERROR_ESTIMATE_HPP
#define DARKSUN_STIFF_RADAU_ERROR_ESTIMATE_HPP

#include "darksun/stiff/common.hpp"
#include "darksun/stiff/decomp.hpp"
#include "darksun/stiff/decsol.hpp"
#include "darksun/stiff/vector_matrix.hpp"
#include <cmath>

namespace darksun::stiff {
template <int N>
void estrad(const Matrix<double, N, N> &fjac, const double h, const double dd1,
            const double dd2, const double dd3, OdeFun<N> &fcn, size_t &nfcn,
            const Vector<double, N> &y0, const Vector<double, N> &y,
            const double x, const Matrix<double, N, N> &e1,
            const Vector<double, N> &z1, const Vector<double, N> &z2,
            const Vector<double, N> &z3, Vector<double, 4 * N> &cont,
            Vector<double, N> &f1, Vector<double, N> &f2,
            const Vector<int, N> &ip1, Vector<double, N> &scal, double &err,
            const bool first, const bool reject) {

  const double hee1 = dd1 / h;
  const double hee2 = dd2 / h;
  const double hee3 = dd3 / h;

  Vector<double, N> contseg{};

  for (int i = 0; i < N; ++i) {
    f2(i) = hee1 * z1(i) + hee2 * z2(i) + hee3 * z3(i);
    cont(i) = f2(i) + y0(i);
  }
  contseg = cont.segment(0, N);
  sol(e1, contseg, ip1);
  cont.segment(0, N) = contseg;

  err = 0.0;
  for (int i = 0; i < N; ++i) {
    /* Computing 2nd power */
    const double d1 = cont(i) / scal(i);
    err += d1 * d1;
  }
  /* Computing MAX */
  double d1 = sqrt(err / double(N));
  err = std::max(d1, 1e-10);

  if (err < 1.0) {
    return;
  }
  if (first || reject) {
    for (int i = 0; i < N; ++i) {
      cont(i) = y(i) + cont(i);
    }
    fcn(x, cont.segment(0, N), f1);
    ++(nfcn);
    for (int i = 0; i < N; ++i) {
      cont(i) = f1(i) + f2(i);
    }
    contseg = cont.segment(0, N);
    sol(e1, contseg, ip1);
    cont.segment(0, N) = contseg;
    err = 0.0;
    for (int i = 0; i < N; ++i) {
      /* Computing 2nd power */
      d1 = cont(i) / scal(i);
      err += d1 * d1;
    }
    /* Computing MAX */
    d1 = sqrt(err / double(N));
    err = std::max(d1, 1e-10);
  }
}

} // namespace darksun::stiff

#endif // DARKSUN_STIFF_RADAU_ERROR_ESTIMATE_HPP
