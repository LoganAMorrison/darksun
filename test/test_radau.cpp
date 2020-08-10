
#include <darksun/stiff/radau5.hpp>
#include <darksun/stiff/vector_matrix.hpp>
#include <fmt/core.h>
#include <gtest/gtest.h>

using namespace darksun::stiff;

TEST(TestRadau5, TestVanDerPol) {

  const double eps = 1e-6;
  constexpr int N = 2;

  OdeFun<N> vdpol = [eps](double, const Vector<double, N> &y,
                          Vector<double, N> &dy) {
    dy(0) = y(1);
    dy(1) = ((1 - y(0) * y(0)) * y(1) - y(0)) / eps;
  };

  OdeJac<N> jvpol = [eps](double, const Vector<double, N> &y,
                          Matrix<double, N, N> &J) {
    J(0, 0) = 0.0;
    J(0, 1) = 1.0;

    J(1, 0) = (-2.0 * y(0) * y(1) - 1.0) / eps;
    J(1, 1) = ((1.0 - y(0) * y(0))) / eps;
  };

  OdeSolOut<N> solout =
      [](const int nr, double xold, double x, const Vector<double, N> &y,
         const Vector<double, 4 * N> &cont, int &irtrn, const Conra5 &conra5) {
        const double y0 = contr5<N>(1, x, cont, conra5);
        const double y1 = contr5<N>(2, x, cont, conra5);

        fmt::print("[{},{},{}],\n", x, y0, y1);
      };

  //=======================
  //---- Initial Conditions
  //=======================
  double xstart = 0.0;
  const double xend = 2.0;
  Vector<double, N> y{};
  y << 2.0, -0.66;

  //=======================
  //---- Tolerances
  //=======================
  Vector<double, N> rtol{}, atol{};
  rtol << 1e-9, 1e-9;
  atol << 1e-9, 1e-9;

  //=======================
  //---- Parameters
  //=======================
  double h = 1e-6;
  const bool iout = true;
  Radau5Parameters params{};
  int idid;

  radau5<N>(vdpol, xstart, y, xend, h, rtol, atol, jvpol, solout, iout, params,
            idid);
}