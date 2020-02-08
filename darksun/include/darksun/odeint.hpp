//
// Created by Logan Morrison on 2/5/20.
//

#pragma once

#include <Eigen/Dense>
#include <limits>
#include <exception>

namespace darksun {

struct Output {
    int kmax;
    int nvar{};
    int nsave{};
    bool dense;
    int count;
    double x1{}, x2{}, xout{}, dxout{};
    Eigen::VectorXd xsave;
    Eigen::MatrixXd ysave;

    Output() : kmax(-1), dense(false), count(0) {}

    explicit Output(const int n_save) : kmax(500), nsave(n_save), count(0), xsave(kmax) {
        dense = n_save > 0;
    }

    void init(const int neqn, const double xlo, const double xhi) {
        nvar = neqn;
        if (kmax == -1)
            return;
        ysave.resize(nvar, kmax);
        if (dense) {
            x1 = xlo;
            x2 = xhi;
            xout = x1;
            dxout = (x2 - x1) / nsave;
        }
    }

    void resize() {
        int kold = kmax;
        kmax *= 2;

        Eigen::VectorXd tempvec(xsave);
        xsave.resize(kmax);
        for (int k = 0; k < kold; k++)
            xsave(k) = tempvec(k);

        Eigen::MatrixXd tempmat(ysave);
        ysave.resize(nvar, kmax);
        for (int i = 0; i < nvar; i++)
            for (int k = 0; k < kold; k++)
                ysave(i, k) = tempmat(i, k);
    }

    template<class Stepper>
    void save_dense(Stepper &s, const double x_out, const double h) {
        if (count == kmax) resize();

        for (int i = 0; i < nvar; i++) ysave(i, count) = s.dense_out(i, x_out, h);
        xsave[count++] = x_out;
    }

    void save(const double x, Eigen::VectorXd &y) {
        if (kmax <= 0) return;

        if (count == kmax) resize();

        for (int i = 0; i < nvar; i++) ysave(i, count) = y(i);

        xsave[count++] = x;
    }

    template<class Stepper>
    void out(const int nstp, const double x, Eigen::VectorXd &y, Stepper &s, const double h) {
        if (!dense) throw std::runtime_error("dense output not set in Output!");

        if (nstp == -1) {
            save(x, y);
            xout += dxout;
        } else {
            while ((x - xout) * (x2 - x1) > 0.0) {
                save_dense(s, xout, h);
                xout += dxout;
            }
        }
    }
};

template<class Stepper>
struct Odeint {
    static const int MAXSTP = 50000;
    double EPS;
    int nok;
    int nbad;
    int nvar;
    double x1, x2, hmin;
    bool dense;
    Eigen::VectorXd y, dydx;
    Eigen::VectorXd &ystart;
    Output &out;
    typename Stepper::Dtype &derivs;
    Stepper s;
    int nstp{};
    double x, h;

    Odeint(Eigen::VectorXd &y_start, double x_1, double x_2, double atol,
           double rtol, double h1, double h_min, Output &output,
           typename Stepper::Dtype &derivatives);

    void integrate();
};

template<class Stepper>
Odeint<Stepper>::Odeint(
        Eigen::VectorXd &y_start, double x_1, double x_2, double atol,
        double rtol, double h1, double h_min, Output &output,
        typename Stepper::Dtype &derivatives)
        : nvar(y_start.size()), y(nvar), dydx(nvar), ystart(y_start), x(x_1), nok(0), nbad(0),
          x1(x_1), x2(x_2), hmin(h_min), dense(output.dense), out(output), derivs(derivatives),
          s(y, dydx, x, atol, rtol, dense) {

    EPS = std::numeric_limits<double>::epsilon();
    h = (x2 - x1 > 0) ? h1 : -h1;

    for (int i = 0; i < nvar; i++) y[i] = ystart[i];

    out.init(s.neqn, x1, x2);
}

template<class Stepper>
void Odeint<Stepper>::integrate() {
    derivs(x, y, dydx);

    if (dense) out.out(-1, x, y, s, h);
    else out.save(x, y);

    for (nstp = 0; nstp < MAXSTP; nstp++) {
        if ((x + h * 1.0001 - x2) * (x2 - x1) > 0.0) h = x2 - x;

        s.step(h, derivs);

        if (s.hdid == h) ++nok;
        else ++nbad;

        if (dense) out.out(nstp, x, y, s, s.hdid);
        else out.save(x, y);

        if ((x - x2) * (x2 - x1) >= 0.0) {
            for (int i = 0; i < nvar; i++) ystart[i] = y[i];

            if (out.kmax > 0 && abs(out.xsave[out.count - 1] - x2) > 100.0 * abs(x2) * EPS)
                out.save(x, y);
            return;
        }
        if (abs(s.hnext) <= hmin) throw std::runtime_error("Step size too small in Odeint");
        h = s.hnext;
    }
    throw std::runtime_error("Too many steps in routine Odeint");
}

}

