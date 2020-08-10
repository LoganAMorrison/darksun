#ifndef DARKSUN_STIFF_RADAU5_HPP
#define DARKSUN_STIFF_RADAU5_HPP

#include "darksun/stiff/common.hpp"
#include "darksun/stiff/decomp.hpp"
#include "darksun/stiff/radau_error_estimate.hpp"
#include "darksun/stiff/radau_solvers.hpp"
#include <iostream>

namespace darksun::stiff {

struct Conra5 {
  size_t nn, nn2, nn3;
  double xsol, hsol, c2m1, c1m1;
};
template <int N>
using OdeSolOut = std::function<void(
    int nr, double xold, double x, const Vector<double, N> &y,
    const Vector<double, 4 * N> &cont, int &irtrn, const Conra5 &contra)>;

struct Radau5Parameters {
  size_t nmax = 100000;
  double uround = std::numeric_limits<double>::epsilon();
  double safe = 0.9;
  double thet = 0.001;
  double fnewt = 0;
  double quot1 = 1.0;
  double quot2 = 1.2;
  double hmax = 0;
  size_t nit = 7;
  bool startn = false;
  bool pred = true;
  double facl = 5.0;
  double facr = 0.125;
};

struct Radau5Statistics {
  size_t nfcn, njac, nstep, naccpt, nrejct, ndec, nsol;
};
template <int N>
double contr5(const int i, const double x, const Vector<double, 4 * N> &cont,
              const Conra5 &conra5) {

  /* ---------------------------------------------------------- */
  /*     THIS FUNCTION CAN BE USED FOR CONINUOUS OUTPUT. IT PROVIDES AN */
  /*     APPROXIMATION TO THE I-TH COMPONENT OF THE SOLUTION AT X. */
  /*     IT GIVES THE VALUE OF THE COLLOCATION POLYNOMIAL, DEFINED FOR */
  /*     THE LAST SUCCESSFULLY COMPUTED STEP (BY RADAU5). */
  /* ---------------------------------------------------------- */

  double s = (x - conra5.xsol) / conra5.hsol;
  return cont(i - 1) + s * (cont(i + conra5.nn - 1) +
                            (s - conra5.c2m1) *
                                (cont(i + conra5.nn2 - 1) +
                                 (s - conra5.c1m1) * cont(i + conra5.nn3 - 1)));
}
template <int N>
int radcor(OdeFun<N> &fcn, double &x, Vector<double, N> &y, const double xend,
           double &h, Vector<double, N> &rtol, Vector<double, N> &atol,
           OdeJac<N> &jac, OdeSolOut<N> &solout, bool iout, int &idid,
           Radau5Parameters &params, Vector<double, N> &z1,
           Vector<double, N> &z2, Vector<double, N> &z3, Vector<double, N> &y0,
           Vector<double, N> &scal, Vector<double, N> &f1,
           Vector<double, N> &f2, Vector<double, N> &f3,
           Matrix<double, N, N> &fjac, Matrix<double, N, N> &e1,
           Matrix<double, N, N> &e2r, Matrix<double, N, N> &e2i,
           Vector<int, N> &ip1, Vector<int, N> &ip2,
           Vector<double, 4 * N> &cont, Radau5Statistics &stats) {

  /* ---------------------------------------------------------- */
  /*     CORE INTEGRATOR FOR RADAU5 */
  /*     PARAMETERS SAME AS IN RADAU5 WITH WORKSPACE ADDED */
  /* ---------------------------------------------------------- */
  /*         DECLARATIONS */
  /* ---------------------------------------------------------- */
  /* *** *** *** *** *** *** *** */
  /*  INITIALISATIONS */
  /* *** *** *** *** *** *** *** */

  Conra5 conra5{};

  /* System generated locals */
  int i1, i_2, i_3, i_4;
  double d1, d_2, d_3, d_4;

  /* Local variables */
  int i, j, k, l;
  double a1, a2, a3;
  int j1;
  size_t n2, n3;
  double ak;
  int md;
  int mm;
  double qt, ak1, ak2, ak3, f1i, f2i, f3i, c1q, c2q, c3q, z1i, z2i, z3i, fac;
  size_t lrc;
  int ier;
  double xph, thq, err, fac1, hacc;
  int lbeg;
  int lend;
  double delt, hnew;
  int newt;
  double dyno, dyth, quot, hhfac, betan, alphn, denom, theta, ysafe;
  int irtrn, nsolu;
  size_t nrsol;
  double qnewt, xosol, acont3;
  bool caljac;

  bool calhes;
  double erracc;
  int mujacj;

  double facgus;
  int mujacp;

  double dynold;

  double thqold;

  /* Function Body */
  conra5.nn = N;
  conra5.nn2 = N * 2;
  conra5.nn3 = N * 3;
  lrc = N << 2;

  /* ---------- CONSTANTS --------- */
  const double sq6 = sqrt(6.);
  const double c1 = (4. - sq6) / 10.;
  const double c2 = (sq6 + 4.) / 10.;
  conra5.c1m1 = c1 - 1.;
  conra5.c2m1 = c2 - 1.;
  const double c1mc2 = c1 - c2;
  const double dd1 = -(sq6 * 7. + 13.) / 3.;
  const double dd2 = (sq6 * 7. - 13.) / 3.;
  const double dd3 = -1.0 / 3.0;
  double u1 = (pow(81.0, (1.0 / 3.0)) + 6. - pow(9.0, (1.0 / 3.0))) / 30.;
  double alph = (12. - pow(81.0, (1.0 / 3.0)) + pow(9.0, (1.0 / 3.0))) / 60.;
  double beta =
      (pow(81.0, (1.0 / 3.0)) + pow(9.0, (1.0 / 3.0))) * sqrt(3.) / 60.;
  /* Computing 2nd power */
  const double cno = alph * alph + beta * beta;
  u1 = 1. / u1;
  alph /= cno;
  beta /= cno;
  const double t11 = .091232394870892942792;
  const double t12 = -.14125529502095420843;
  const double t13 = -.030029194105147424492;
  const double t21 = .24171793270710701896;
  const double t22 = .20412935229379993199;
  const double t23 = .38294211275726193779;
  const double t31 = .96604818261509293619;
  const double ti11 = 4.325579890063155351;
  const double ti12 = .33919925181580986954;
  const double ti13 = .54177053993587487119;
  const double ti21 = -4.1787185915519047273;
  const double ti22 = -.32768282076106238708;
  const double ti23 = .47662355450055045196;
  const double ti31 = -.50287263494578687595;
  const double ti32 = 2.5719269498556054292;
  const double ti33 = -.59603920482822492497;

  double posneg = std::copysign(1.0, xend - x);
  /* Computing MIN */
  double hmaxn = std::min(std::abs(params.hmax), std::abs(xend - x));
  if (abs(h) <= params.uround * 10.) {
    h = 1e-6;
  }
  /* Computing MIN */
  h = std::min(std::abs(h), hmaxn);
  h = std::copysign(h, posneg);
  double hold = h;
  bool reject = false;
  bool first = true;
  bool last = false;
  if ((x + h * 1.0001 - xend) * posneg >= 0.) {
    h = xend - x;
    last = true;
  }
  double hopt = h;
  double faccon = 1.;
  double cfac = params.safe * double((params.nit << size_t(1)) + 1);
  int nsing = 0;
  double xold = x;
  if (iout) {
    irtrn = 1;
    nrsol = 1;
    xosol = xold;
    conra5.xsol = x;
    for (i = 1; i <= N; ++i) {
      cont(i - 1) = y(i - 1);
    }
    nsolu = N;
    conra5.hsol = hold;
    solout(nrsol, xosol, conra5.xsol, y, cont, irtrn, conra5);
    if (irtrn < 0) {
      goto L179;
    }
  }
  n2 = N * 2;
  n3 = N * 3;
  for (i = 1; i <= N; ++i) {
    scal(i - 1) = atol(i - 1) + rtol(i - 1) * (d1 = y(i - 1), abs(d1));
  }
  hhfac = h;
  fcn(x, y, y0);
  ++stats.nfcn;
/* --- BASIC INTEGRATION STEP */
L10:
  /* *** *** *** ** *** *** *** */
  /*  COMPUTATION OF THE JACOBIAN */
  /* *** *** *** *** *** *** *** */
  ++stats.njac;
  jac(x, y, fjac);
  caljac = true;
  calhes = true;
L20:
  /* --- COMPUTE THE MATRICES E1 AND E2 AND THEIR DECOMPOSITIONS */
  fac1 = u1 / h;
  alphn = alph / h;
  betan = beta / h;
  ier = decomr(fjac, fac1, e1, ip1);
  if (ier != 0) {
    goto L78;
  }
  ier = decomc(fjac, alphn, betan, e2r, e2i, ip2);
  if (ier != 0) {
    goto L78;
  }
  ++stats.ndec;
L30:
  ++stats.nstep;
  if (stats.nstep > params.nmax) {
    goto L178;
  }
  if (abs(h) * 0.1 <= abs(x) * params.uround) {
    goto L177;
  }
  xph = x + h;
  /* *** *** *** *** *** *** *** */
  /*  STARTING VALUES FOR NEWTON ITERATION */
  /* *** *** *** *** *** *** *** */
  if (first || params.startn) {
    for (i = 1; i <= N; ++i) {
      z1(i - 1) = 0.;
      z2(i - 1) = 0.;
      z3(i - 1) = 0.;
      f1(i - 1) = 0.;
      f2(i - 1) = 0.;
      f3(i - 1) = 0.;
    }
  } else {
    c3q = h / hold;
    c1q = c1 * c3q;
    c2q = c2 * c3q;
    for (i = 1; i <= N; ++i) {
      ak1 = cont(i + N - 1);
      ak2 = cont(i + n2 - 1);
      ak3 = cont(i + n3 - 1);
      z1i =
          c1q * (ak1 + (c1q - conra5.c2m1) * (ak2 + (c1q - conra5.c1m1) * ak3));
      z2i =
          c2q * (ak1 + (c2q - conra5.c2m1) * (ak2 + (c2q - conra5.c1m1) * ak3));
      z3i =
          c3q * (ak1 + (c3q - conra5.c2m1) * (ak2 + (c3q - conra5.c1m1) * ak3));
      z1(i - 1) = z1i;
      z2(i - 1) = z2i;
      z3(i - 1) = z3i;
      f1(i - 1) = ti11 * z1i + ti12 * z2i + ti13 * z3i;
      f2(i - 1) = ti21 * z1i + ti22 * z2i + ti23 * z3i;
      f3(i - 1) = ti31 * z1i + ti32 * z2i + ti33 * z3i;
    }
  }
  /* *** *** *** *** *** *** *** */
  /*  LOOP FOR THE SIMPLIFIED NEWTON ITERATION */
  /* *** *** *** *** *** *** *** */
  newt = 0;
  d1 = std::max(faccon, params.uround);
  faccon = pow(d1, 0.8);
  theta = abs(params.thet);
L40:
  if (newt >= params.nit) {
    goto L78;
  }
  /* ---     COMPUTE THE RIGHT-HAND SIDE */
  for (i = 1; i <= N; ++i) {
    cont(i - 1) = y(i - 1) + z1(i - 1);
  }
  d1 = x + c1 * h;
  fcn(d1, cont.segment(0, N), z1);
  for (i = 1; i <= N; ++i) {
    cont(i - 1) = y(i - 1) + z2(i - 1);
  }
  d1 = x + c2 * h;
  fcn(d1, cont.segment(0, N), z2);
  for (i = 1; i <= N; ++i) {
    cont(i - 1) = y(i - 1) + z3(i - 1);
  }
  fcn(xph, cont.segment(0, N), z3);
  stats.nfcn += 3;
  /* ---     SOLVE THE LINEAR SYSTEMS */
  for (i = 1; i <= N; ++i) {
    a1 = z1(i - 1);
    a2 = z2(i - 1);
    a3 = z3(i - 1);
    z1(i - 1) = ti11 * a1 + ti12 * a2 + ti13 * a3;
    z2(i - 1) = ti21 * a1 + ti22 * a2 + ti23 * a3;
    z3(i - 1) = ti31 * a1 + ti32 * a2 + ti33 * a3;
  }

  slvrad(fjac, fac1, alphn, betan, e1, e2r, e2i, z1, z2, z3, f1, f2, f3, ip1,
         ip2);

  ++stats.nsol;
  ++newt;
  dyno = 0.;
  for (i = 1; i <= N; ++i) {
    denom = scal(i - 1);
    /* Computing 2nd power */
    d1 = z1(i - 1) / denom;
    /* Computing 2nd power */
    d_2 = z2(i - 1) / denom;
    /* Computing 2nd power */
    d_3 = z3(i - 1) / denom;
    dyno = dyno + d1 * d1 + d_2 * d_2 + d_3 * d_3;
  }
  dyno = sqrt(dyno / n3);
  /* ---     BAD CONVERGENCE OR NUMBER OF ITERATIONS TO LARGE */
  if (newt > 1 && newt < params.nit) {
    thq = dyno / dynold;
    if (newt == 2) {
      theta = thq;
    } else {
      theta = sqrt(thq * thqold);
    }
    thqold = thq;
    if (theta < .99) {
      faccon = theta / (1. - theta);
      i1 = int(params.nit) - 1 - newt;
      dyth = faccon * dyno * pow(theta, i1) / params.fnewt;
      if (dyth >= 1.) {
        /* Computing MAX */
        d1 = 1e-4, d_2 = std::min(20., dyth);
        qnewt = std::max(d1, d_2);
        d1 = -1. / (params.nit + 4.0 - 1 - newt);
        hhfac = pow(qnewt, d1) * 0.8;
        h = hhfac * h;
        reject = true;
        last = false;
        if (caljac) {
          goto L20;
        }
        goto L10;
      }
    } else {
      goto L78;
    }
  }
  dynold = std::max(dyno, params.uround);
  for (i = 1; i <= N; ++i) {
    f1i = f1(i - 1) + z1(i - 1);
    f2i = f2(i - 1) + z2(i - 1);
    f3i = f3(i - 1) + z3(i - 1);
    f1(i - 1) = f1i;
    f2(i - 1) = f2i;
    f3(i - 1) = f3i;
    z1(i - 1) = t11 * f1i + t12 * f2i + t13 * f3i;
    z2(i - 1) = t21 * f1i + t22 * f2i + t23 * f3i;
    z3(i - 1) = t31 * f1i + f2i;
  }
  if (faccon * dyno > params.fnewt) {
    goto L40;
  }
  /* --- ERROR ESTIMATION */
  estrad(fjac, h, dd1, dd2, dd3, fcn, stats.nfcn, y0, y, x, e1, z1, z2, z3,
         cont, f1, f2, ip1, scal, err, first, reject);

  /* --- COMPUTATION OF HNEW */
  /* --- WE REQUIRE .2<=HNEW/H<=8. */
  /* Computing MIN */
  d1 = params.safe, d_2 = cfac / double(newt + (params.nit << size_t(1)));
  fac = std::min(d1, d_2);
  /* Computing MAX */
  /* Computing MIN */
  d_3 = params.facl, d_4 = pow(err, 0.25) / fac;
  d1 = params.facr, d_2 = std::min(d_3, d_4);
  quot = std::max(d1, d_2);
  hnew = h / quot;
  /* *** *** *** *** *** *** *** */
  /*  IS THE ERROR SMALL ENOUGH ? */
  /* *** *** *** *** *** *** *** */
  if (err < 1.) {
    /* --- STEP IS ACCEPTED */
    first = false;
    ++stats.naccpt;
    if (params.pred) {
      /*       --- PREDICTIVE CONTROLLER OF GUSTAFSSON */
      if (stats.naccpt > 1) {
        /* Computing 2nd power */
        d_2 = err;
        d1 = d_2 * d_2 / erracc;
        facgus = hacc / h * pow(d1, 0.25) / params.safe;
        /* Computing MAX */
        d1 = params.facr, d_2 = std::min(params.facl, facgus);
        facgus = std::max(d1, d_2);
        quot = std::max(quot, facgus);
        hnew = h / quot;
      }
      hacc = h;
      erracc = std::max(.01, err);
    }
    xold = x;
    hold = h;
    x = xph;
    for (i = 1; i <= N; ++i) {
      y(i - 1) += z3(i - 1);
      z2i = z2(i - 1);
      z1i = z1(i - 1);
      cont(i + N - 1) = (z2i - z3(i - 1)) / conra5.c2m1;
      ak = (z1i - z2i) / c1mc2;
      acont3 = z1i / c1;
      acont3 = (ak - acont3) / c2;
      cont(i + n2 - 1) = (ak - cont(i + N - 1)) / conra5.c1m1;
      cont(i + n3 - 1) = cont(i + n2 - 1) - acont3;
    }
    for (i = 1; i <= N; ++i) {
      scal(i - 1) = atol(i - 1) + rtol(i - 1) * (d1 = y(i - 1), abs(d1));
    }
    if (iout) {
      nrsol = stats.naccpt + 1;
      conra5.xsol = x;
      xosol = xold;
      for (i = 1; i <= N; ++i) {
        cont(i - 1) = y(i - 1);
      }
      nsolu = N;
      conra5.hsol = hold;
      solout(nrsol, xosol, conra5.xsol, y, cont, irtrn, conra5);
      if (irtrn < 0) {
        goto L179;
      }
    }
    caljac = false;
    if (last) {
      h = hopt;
      idid = 1;
      return 0;
    }
    fcn(x, y, y0);
    ++stats.nfcn;
    /* Computing MIN */
    d1 = abs(hnew);
    hnew = posneg * std::min(d1, hmaxn);
    hopt = hnew;
    hopt = std::min(h, hnew);
    if (reject) {
      /* Computing MIN */
      d1 = abs(hnew), d_2 = abs(h);
      hnew = posneg * std::min(d1, d_2);
    }
    reject = false;
    if ((x + hnew / params.quot1 - xend) * posneg >= 0.) {
      h = xend - x;
      last = true;
    } else {
      qt = hnew / h;
      hhfac = h;
      if (theta <= params.thet && qt >= params.quot1 && qt <= params.quot2) {
        goto L30;
      }
      h = hnew;
    }
    hhfac = h;
    if (theta <= params.thet) {
      goto L20;
    }
    goto L10;
  } else {
    /* --- STEP IS REJECTED */
    reject = true;
    last = false;
    if (first) {
      h *= .1;
      hhfac = .1;
    } else {
      hhfac = hnew / h;
      h = hnew;
    }
    if (stats.naccpt >= 1) {
      ++stats.nrejct;
    }
    if (caljac) {
      goto L20;
    }
    goto L10;
  }
/* --- UNEXPECTED STEP-REJECTION */
L78:
  if (ier != 0) {
    ++nsing;
    if (nsing >= 5) {
      goto L176;
    }
  }
  h *= .5;
  hhfac = .5;
  reject = true;
  last = false;
  if (caljac) {
    goto L20;
  }
  goto L10;
/* --- FAIL EXIT */
L176:
  std::cerr << "Encountered repeatedly singular matrix, ier = " << ier << "\n";
  idid = -4;
  return 0;
L177:
  std::cerr << "Step size too small: h = " << h << "\n";
  idid = -3;
  return 0;
L178:
  std::cerr << " More than nmax = " << params.nmax << " steps are needed.\n";
  idid = -2;
  return 0;
/* --- EXIT CAUSED BY SOLOUT */
L179:
  idid = 2;
  return 0;
}
template <int N>
int radau5(OdeFun<N> &fcn, double &x, Vector<double, N> &y, double xend,
           double h, Vector<double, N> &rtol, Vector<double, N> &atol,
           OdeJac<N> &jac, OdeSolOut<N> &solout, bool iout,
           Radau5Parameters &params, int &idid) {

  Radau5Statistics stats{};

  if (params.uround <= 1e-19 || params.uround >= 1.) {
    std::cerr << "Coefficients have 20 digits: uround = " << params.uround
              << "\n";
  }
  /* -------- CHECK AND CHANGE THE TOLERANCES */
  double expm = .66666666666666663;
  for (int i = 1; i <= N; ++i) {
    if (atol(i - 1) <= 0. || rtol(i - 1) <= params.uround * 10.) {
      throw std::runtime_error("Tolerances(" + std::to_string(i) +
                               ") are too small.");
    } else {
      double quot = atol(i - 1) / rtol(i - 1);
      rtol(i - 1) = pow(rtol(i - 1), expm) * 0.1;
      atol(i - 1) = rtol(i - 1) * quot;
    }
  }
  /* -------- NMAX , THE MAXIMAL NUMBER OF STEPS ----- */
  if (params.nmax == 0) {
    params.nmax = 100000;
  } else {
    if (params.nmax <= 0) {
      throw std::runtime_error("Invalid value: nmax = " +
                               std::to_string(params.nmax));
    }
  }
  /* -------- NIT    MAXIMAL NUMBER OF NEWTON ITERATIONS */
  if (params.nit <= 0) {
    throw std::runtime_error("Invalid value: nit = " +
                             std::to_string(params.nit));
  }
  /* --------- SAFE     SAFETY FACTOR IN STEP SIZE PREDICTION */
  if (params.safe <= .001 || params.safe >= 1.) {
    throw std::runtime_error("Invalid value: safe = " +
                             std::to_string(params.safe));
  }
  /* ------ THET     DECIDES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED; */
  if (params.thet >= 1.) {
    throw std::runtime_error("Invalid value: thet = " +
                             std::to_string(params.thet));
  }
  /* --- FNEWT   STOPPING CRITERION FOR NEWTON'S METHOD, USUALLY CHOSEN <1. */
  double tolst = rtol(1 - 1);
  if (params.fnewt == 0.) {
    /* Computing MAX */
    /* Computing MIN */
    params.fnewt =
        std::max(params.uround * 10 / tolst, std::min(0.03, pow(tolst, 0.5)));
  } else {
    if (params.fnewt <= params.uround / tolst) {
      throw std::runtime_error("Invalid value: fnewt = " +
                               std::to_string(params.fnewt));
    }
  }
  /* --- QUOT1 AND QUOT2: IF QUOT1 < HNEW/HOLD < QUOT2, STEP SIZE = CONST. */
  if (params.quot1 > 1. || params.quot2 < 1.) {
    throw std::runtime_error(
        "Invalid value(s): quot1, quot2 = " + std::to_string(params.quot1) +
        ", " + std::to_string(params.quot2));
  }
  /* -------- MAXIMAL STEP SIZE */
  if (params.hmax == 0) {
    params.hmax = xend - x;
  }
  /* -------  FACL,FACR     PARAMETERS FOR STEP SIZE SELECTION */
  if (params.facl < 1.0 || params.facr > 1.0) {
    throw std::runtime_error(
        "Invalid value(s): facl, facr = " + std::to_string(params.facl) + ", " +
        std::to_string(params.facr));
  }

  //===============================================================================
  //---- Array creation
  //===============================================================================
  Vector<double, N> z1{};
  Vector<double, N> z2{};
  Vector<double, N> z3{};
  Vector<double, N> f1{};
  Vector<double, N> f2{};
  Vector<double, N> f3{};
  Vector<double, N> y0{};
  Vector<double, N> scal{};
  Vector<double, 4 * N> cont{};

  Vector<int, N> ip1{};
  Vector<int, N> ip2{};

  Matrix<double, N, N> fjac{};
  Matrix<double, N, N> e1{};
  Matrix<double, N, N> e2r{};
  Matrix<double, N, N> e2i{};

  radcor(fcn, x, y, xend, h, rtol, atol, jac, solout, iout, idid, params, z1,
         z2, z3, y0, scal, f1, f2, f3, fjac, e1, e2r, e2i, ip1, ip2, cont,
         stats);
  /* -------- RESTORE TOLERANCES */
  expm = 1. / expm;
  for (int i = 1; i < N; ++i) {
    double quot = atol(i - 1) / rtol(i - 1);
    double d1 = rtol(i - 1) * 10.;
    rtol(i - 1) = pow(d1, expm);
    atol(i - 1) = rtol(i - 1) * quot;
  }
  return 0;
}

} // namespace darksun::stiff

#endif // DARKSUN_STIFF_RADAU5_HPP
