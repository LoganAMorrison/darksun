
#include <fmt/core.h>
#include <gtest/gtest.h>
#include <stiff/stiff.hpp>

using namespace stiff;

template <int N, int M> struct Solution {
  std::array<double, M> x;
  std::array<std::array<double, N>, M> y;
};

TEST(TestRadau, TestVanDerPol) {

  int nd = 2;     // dimension of the system
  int ns = 7;     // maximum number of allowed stages
  int ijac = 1;   // analytic Jacobian
  int mljac = nd; // full Jacaobian
  int mujac = 0;  // full Jacaobian
  int imas = 0;   // no mass matrix
  int mlmas = 0;  // no mass matrix
  int mumas = 0;  // no mass matrix

  // Length of `work`
  int lwork = (ns + 1) * nd * nd + (3 * ns + 3) * nd + 20;
  // Length of `iwork`
  int liwork = (2 + (ns - 1) / 2) * nd + 20;
  std::vector<double> work(lwork); // Workspace of doubles needed for radau
  std::vector<int> iwork(liwork);  // Workspace of ints needed for radau
  for (int i = 0; i < 20; i++) {   // Set all radau params to defaults
    iwork[i] = 0;
    work[i] = 0.0;
  }

  int idid;                 // Flag for radau return code
  double rtol = 1.0e-9;     // Relative tolerance for radau
  double atol = 1.0 * rtol; // Absolute tolerance for radau
  int itol = 0;             // Flag telling radau we are using scalar tols
  double h = 1.0e-6;        // Initial step size
  int iout = 1;             // Have radau call `solout`

  double xstart = 0.0;
  double xend = 2.0;
  double y[2];
  y[0] = 2.0;
  y[1] = -0.66;

  constexpr int num_xs = 41;
  Solution<2, 41> sol{};
  double eps = 1e-6;
  double dxout = (xend - xstart) / double(num_xs - 1);

  auto vdpol = [eps](int *, double *, double *y, double *dy) {
    dy[0] = y[1];
    dy[1] = ((1 - y[0] * y[0]) * y[1] - y[0]) / eps;
  };

  auto jvpol = [eps](int *, double *, double *y, double *dfy, int *) {
    dfy[0] = 0.0;
    dfy[1] = (-2.0 * y[0] * y[1] - 1.0) / eps;
    dfy[2] = 1.0;
    dfy[3] = ((1 - y[0] * y[0])) / eps;
  };
  auto mas = [](int *, double *, int *) {};
  auto solout = [dxout, &sol](int *nr, double *xold, double *x, double *,
                              double *cont, int *lrc, int *, int *,
                              RadauWeight &w) {
    static double d = 0;
    static size_t i = 0;
    double yd[2];
    if (*nr == 1) {
      d = *xold;
    }

    int idx1 = 1, idx2 = 2;
    yd[0] = contra(&idx1, x, cont, lrc, w);
    yd[1] = contra(&idx2, x, cont, lrc, w);
    fmt::print("[{},{},{}],\n", *x, yd[0], yd[1]);
    while ((*xold <= d) && (*x >= d)) {
      yd[0] = contra(&idx1, &d, cont, lrc, w);
      yd[1] = contra(&idx2, &d, cont, lrc, w);
      sol.x[i] = d;
      sol.y[i][0] = yd[0];
      sol.y[i][1] = yd[1];

      d += dxout;
      i += 1;
    }
  };

  radau(&nd, vdpol, &xstart, y, &xend, &h, &rtol, &atol, &itol, jvpol, &ijac,
        &mljac, &mujac, mas, &imas, &mlmas, &mumas, solout, &iout, work.data(),
        &lwork, iwork.data(), &liwork, &idid);
  // Save final solution
  sol.x[(num_xs - 1)] = xend;
  sol.y[(num_xs - 1)][0] = y[0];
  sol.y[(num_xs - 1)][1] = y[1];

  //=========================================================================
  //---- Mathematica Solution -----------------------------------------------
  //=========================================================================

  std::vector<std::vector<double>> mma{
      {0., 2., -0.66},
      {0.05, 1.9661892484754215, -0.6860638715110897},
      {0.1, 1.9313610478721095, -0.7074178507553899},
      {0.15000000000000002, 1.8954088384677206, -0.7310918015640505},
      {0.2, 1.8582056276363659, -0.7575455973910044},
      {0.25, 1.8195979561669557, -0.7873851758826808},
      {0.30000000000000004, 1.7793974167757798, -0.8214160097192423},
      {0.35000000000000003, 1.7373683162108742, -0.8607441151327514},
      {0.4, 1.6932090256578511, -0.9069350923615699},
      {0.45, 1.6465226506589787, -0.9622955673428706},
      {0.5, 1.5967684573030676, -1.0303919311654355},
      {0.55, 1.543175413976293, -1.117116105141874},
      {0.6000000000000001, 1.4845750932363537, -1.2330713461728815},
      {0.65, 1.4190350204854123, -1.399907588330944},
      {0.7000000000000001, 1.342890898941202, -1.6715916588912154},
      {0.75, 1.2472019295857966, -2.2451024277919465},
      {0.8, 1.0839201216017382, -6.195465640457782},
      {0.8500000000000001, -1.971118662159346, 0.6831570250070442},
      {0.9, -1.936443129338714, 0.7042090923029367},
      {0.9500000000000001, -1.9006600748947884, 0.7275230893444639},
      {1., -1.8636457734095244, 0.7535434289873758},
      {1.05, -1.8252510347857684, 0.7828513484575084},
      {1.1, -1.7852931284223588, 0.8162189064954589},
      {1.1500000000000001, -1.7435441383035513, 0.8547008616567286},
      {1.2000000000000002, -1.6997134197726316, 0.8997826986496342},
      {1.25, -1.6534202018998279, 0.9536402719476197},
      {1.3, -1.6041486131673408, 1.0196113464070558},
      {1.35, -1.5511684254192748, 1.1031511437264612},
      {1.4000000000000001, -1.4933841570646318, 1.2139377356322454},
      {1.4500000000000002, -1.4290107372008636, 1.3713136172955742},
      {1.5, -1.3547447410455409, 1.6217935234704024},
      {1.55, -1.2629136764948916, 2.122693038582543},
      {1.6, -1.1208100241551344, 4.373902091600681},
      {1.6500000000000001, 1.976027273175076, -0.6802900790098091},
      {1.7000000000000002, 1.9415022557335877, -0.7010472435573056},
      {1.75, 1.9058857983366582, -0.7240101812979894},
      {1.8, 1.8690573261425005, -0.7496089960511384},
      {1.85, 1.8308717492878046, -0.7784013725146435},
      {1.9000000000000001, 1.791151777035892, -0.8111271853438429},
      {1.9500000000000002, 1.7496769138953059, -0.8487931530970122},
      {2., 1.706166942293804, -0.8928105392359631}};

  for (size_t i = 0; i < num_xs; i++) {
    double x = sol.x[i];
    double y0 = sol.y[i][0];
    double y1 = sol.y[i][1];
    fmt::print("x = {}, y = ({}, {})\n", x, y0, y1);
    if (mma[i][0] != 0.0) {
      ASSERT_LE(abs((mma[i][0] - x) / mma[i][0]), 1e-5);
    }
    ASSERT_LE(abs((mma[i][1] - y0) / mma[i][1]), 1e-4);
    ASSERT_LE(abs((mma[i][2] - y1) / mma[i][2]), 1e-4);
  }
}