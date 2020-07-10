use super::linalg::*;
use ndarray::prelude::*;

pub struct Radau5 {
    pub n: usize,       // Dimension of the system
    pub y: Array1<f64>, // state of the system
    pub x: f64,         // independent variable
    pub xend: f64,      // final value of independent variable
    pub dx: f64,        // time step for intermediate output
    pub h: f64,         // step size
    // Generic Options
    pub(crate) rtoler: Vec<f64>, // vector of relative tolerances
    pub(crate) atoler: Vec<f64>, // vector of absolute tolerances
    pub(crate) hmax: f64,        // max step size,
    pub(crate) nmax: usize,      // max number of steps
    pub(crate) safe: f64,        // safety factor for step size prediction
    pub(crate) facl: f64,        // parameter for step size selection
    pub(crate) facr: f64,        // parameter for step size selection
    pub(crate) xold: f64,        // stores past value of x
    pub(crate) hold: f64,        // stores past value of h
    pub(crate) xd: f64,          // x at discrete points specified by dx interval
    // Radau options
    pub(crate) nit: usize,   // maximum number of newton iterations
    pub(crate) startn: bool, // switch for starting values of newton iterations
    pub(crate) pred: bool,   // step size control
    pub(crate) fnewt: f64,   // stopping criterion for Newton's method
    pub(crate) quot1: f64,   // if quot1 < hnew/hold < quot2, step size = const
    pub(crate) quot2: f64,   // if quot1 < hnew/hold < quot2, step size = const
    pub(crate) thet: f64,    // decides whether the Jacobian should be recomputed
    // Cache
    pub(crate) naccpt: usize,
    pub(crate) nrejct: usize,
    pub(crate) nstep: usize,
    pub(crate) ndec: usize,
    pub(crate) nfcn: usize,
    pub(crate) nsol: usize,
    pub(crate) fac1: f64,
    pub(crate) alphn: f64,
    pub(crate) betan: f64,
    pub(crate) err: f64,
    pub(crate) caljac: bool,
    pub(crate) first: bool,
    pub(crate) reject: bool,
    pub(crate) z1: Array1<f64>,
    pub(crate) z2: Array1<f64>,
    pub(crate) z3: Array1<f64>,
    pub(crate) f1: Array1<f64>,
    pub(crate) f2: Array1<f64>,
    pub(crate) f3: Array1<f64>,
    pub(crate) y0: Array1<f64>,
    pub(crate) scal: Array1<f64>,
    pub(crate) cont: Array1<f64>,
    pub(crate) ip1: Array1<i32>,
    pub(crate) ip2: Array1<i32>,
    pub(crate) e1: Array2<f64>,
    pub(crate) e2r: Array2<f64>,
    pub(crate) e2i: Array2<f64>,
    pub(crate) fjac: Array2<f64>,
}

impl Radau5 {
    pub fn builder(y: Array1<f64>, xspan: (f64, f64)) -> Radau5Builder {
        Radau5Builder {
            n: y.len(),
            y,
            x: xspan.0,
            xend: xspan.1,
            h: None,
            dx: None,
            rtoler: None,
            atoler: None,
            hmax: None,
            nmax: None,
            safe: None,
            facl: None,
            facr: None,
            xold: None,
            hold: None,
            xd: None,
            nit: None,
            startn: None,
            pred: None,
            fnewt: None,
            quot1: None,
            quot2: None,
            thet: None,
        }
    }
}

pub struct Radau5Builder {
    pub n: usize,       // Dimension of the system
    pub y: Array1<f64>, // state of the system
    pub x: f64,         // independent variable
    pub xend: f64,      // final value of independent variable
    // Generic Options
    pub h: Option<f64>,           // step size
    pub dx: Option<f64>,          // time step for intermediate output
    pub rtoler: Option<Vec<f64>>, // vector of relative tolerances
    pub atoler: Option<Vec<f64>>, // vector of absolute tolerances
    pub hmax: Option<f64>,        // max step size,
    pub nmax: Option<usize>,      // max number of steps
    pub safe: Option<f64>,        // safety factor for step size prediction
    pub facl: Option<f64>,        // parameter for step size selection
    pub facr: Option<f64>,        // parameter for step size selection
    pub xold: Option<f64>,        // stores past value of x
    pub hold: Option<f64>,        // stores past value of h
    pub xd: Option<f64>,          // x at discrete points specified by dx interval
    // Radau options
    pub nit: Option<usize>,   // maximum number of newton iterations
    pub startn: Option<bool>, // switch for starting values of newton iterations
    pub pred: Option<bool>,   // step size control
    pub fnewt: Option<f64>,   // stopping criterion for Newton's method
    pub quot1: Option<f64>,   // if quot1 < hnew/hold < quot2, step size = const
    pub quot2: Option<f64>,   // if quot1 < hnew/hold < quot2, step size = const
    pub thet: Option<f64>,    // decides whether the Jacobian should be recomputed
}

impl Radau5Builder {
    pub fn h(mut self, h: f64) -> Self {
        self.h = Some(h);
        self
    }
    pub fn dx(mut self, dx: f64) -> Self {
        self.dx = Some(dx);
        self
    }
    pub fn rtoler(mut self, rtoler: Vec<f64>) -> Self {
        self.rtoler = Some(rtoler);
        self
    }
    pub fn atoler(mut self, atoler: Vec<f64>) -> Self {
        self.atoler = Some(atoler);
        self
    }
    pub fn hmax(mut self, hmax: f64) -> Self {
        self.hmax = Some(hmax);
        self
    }
    pub fn nmax(mut self, nmax: usize) -> Self {
        self.nmax = Some(nmax);
        self
    }
    pub fn safe(mut self, safe: f64) -> Self {
        self.safe = Some(safe);
        self
    }
    pub fn facl(mut self, facl: f64) -> Self {
        self.facl = Some(facl);
        self
    }
    pub fn facr(mut self, facr: f64) -> Self {
        self.facr = Some(facr);
        self
    }
    pub fn nit(mut self, nit: usize) -> Self {
        self.nit = Some(nit);
        self
    }
    pub fn startn(mut self, startn: bool) -> Self {
        self.startn = Some(startn);
        self
    }
    pub fn pred(mut self, pred: bool) -> Self {
        self.pred = Some(pred);
        self
    }
    pub fn fnewt(mut self, fnewt: f64) -> Self {
        self.fnewt = Some(fnewt);
        self
    }
    pub fn quot1(mut self, quot1: f64) -> Self {
        self.quot1 = Some(quot1);
        self
    }
    pub fn quot2(mut self, quot2: f64) -> Self {
        self.quot2 = Some(quot2);
        self
    }
    pub fn thet(mut self, thet: f64) -> Self {
        self.thet = Some(thet);
        self
    }
    pub fn build(self) -> Result<Radau5, String> {
        let n = self.n;

        // Maximum stepsize
        let hmax = self.hmax.unwrap_or(self.xend - self.x);

        // Max number of steps
        let nmax = self.nmax.unwrap_or(100000);
        if nmax == 0 {
            return Err(format!("Invalid input: nmax = {}", nmax));
        }

        // Safety factor for step size prediction
        let safe = self.safe.unwrap_or(0.9);
        if safe <= 0.001 || safe >= 1.0 {
            return Err(format!("Invalid input: safe = {}", safe));
        }

        // Check and change tolerances:
        let mut rtoler = self.rtoler.unwrap_or(vec![1e-7; self.n]);
        let mut atoler = self.atoler.unwrap_or(vec![1e-7; self.n]);

        for (i, (at, rt)) in atoler.iter_mut().zip(rtoler.iter_mut()).enumerate() {
            if *at <= 0.0 || *rt <= 10.0 * f64::EPSILON {
                return Err(format!(
                    "Tolerances are too small atol[{}]={:?}, rtol[{}]={:?}",
                    i, at, i, rt
                ));
            } else {
                let quot = *at / *rt;
                *rt = 0.1 * (*rt).powf(2.0 / 3.0);
                *at = *rt * quot;
            }
        }

        // Initial step length
        let mut h = self.h.unwrap_or(0.0);
        if h.abs() < 10.0 * f64::EPSILON {
            h = 1e-6;
        }

        // Set facr and facl for step-size selection
        let facl = self.facl.unwrap_or(5.0);
        let facr = self.facr.unwrap_or(1.0 / 8.0);
        if facl < 1.0 || facr > 1.0 {
            return Err(format!("Invalid inputs: facl, facr = {}, {}", facl, facr));
        }

        // Number of newton iterations
        let nit = self.nit.unwrap_or(7);
        if nit <= 0 {
            return Err(format!("Invalid input: nit = {}", nit));
        }

        // Swith for newton iters
        let startn = self.startn.unwrap_or(false);

        // Step size predition
        let pred = self.pred.unwrap_or(true);

        // Stopping criterion for Netwon's method
        let rtolmin = rtoler
            .iter()
            .fold(rtoler[0], |acc, x| if *x < acc { *x } else { acc });
        let fnewt = self
            .fnewt
            .unwrap_or((10.0 * f64::EPSILON / rtolmin).max(0.03f64.min(rtolmin.sqrt())));
        if fnewt <= f64::EPSILON / rtolmin {
            return Err(format!("Invalid input: fnewt = {}", fnewt));
        }

        // Step size controllers
        let quot1 = self.quot1.unwrap_or(1.0);
        let quot2 = self.quot2.unwrap_or(1.2);
        if quot1 > 1.0 || quot2 < 1.0 {
            return Err(format!(
                "Invalid inputs: quot1, quot2 = {}, {}",
                quot1, quot2
            ));
        }

        // thet--decides whether the Jacobian should be recomputed
        let thet = self.thet.unwrap_or(0.001);
        if thet >= 1.0 {
            return Err(format!("Invalid input: thet = {}", thet));
        }

        Ok(Radau5 {
            n: self.n,
            y: self.y,
            x: self.x,
            xend: self.xend,
            dx: self.dx.unwrap_or(hmax),
            h,
            naccpt: 0,
            nrejct: 0,
            nstep: 0,
            ndec: 0,
            nfcn: 0,
            nsol: 0,
            rtoler,
            atoler,
            hmax,
            nmax,
            safe,
            facl,
            facr,
            xold: self.x,
            hold: h,
            xd: self.x,
            nit,
            startn,
            pred,
            fnewt,
            quot1,
            quot2,
            thet,
            fac1: 0.0,
            alphn: 0.0,
            betan: 0.0,
            err: 0.0,
            caljac: true,
            first: true,
            reject: false,
            z1: Array1::zeros(n),
            z2: Array1::zeros(n),
            z3: Array1::zeros(n),
            f1: Array1::zeros(n),
            f2: Array1::zeros(n),
            f3: Array1::zeros(n),
            y0: Array1::zeros(n),
            scal: Array1::zeros(n),
            cont: Array1::zeros(4 * n),
            ip1: Array1::zeros(n),
            ip2: Array1::zeros(n),
            e1: Array2::zeros((n, n)),
            e2r: Array2::zeros((n, n)),
            e2i: Array2::zeros((n, n)),
            fjac: Array2::zeros((n, n)),
        })
    }
}

// Declare all the constants needed for Radau5
impl Radau5 {
    pub(crate) const T11: f64 = 9.1232394870892942792e-02;
    pub(crate) const T12: f64 = -0.14125529502095420843;
    pub(crate) const T13: f64 = -3.0029194105147424492e-02;
    pub(crate) const T21: f64 = 0.24171793270710701896;
    pub(crate) const T22: f64 = 0.20412935229379993199;
    pub(crate) const T23: f64 = 0.38294211275726193779;
    pub(crate) const T31: f64 = 0.96604818261509293619;
    pub(crate) const TI11: f64 = 4.3255798900631553510;
    pub(crate) const TI12: f64 = 0.33919925181580986954;
    pub(crate) const TI13: f64 = 0.54177053993587487119;
    pub(crate) const TI21: f64 = -4.1787185915519047273;
    pub(crate) const TI22: f64 = -0.32768282076106238708;
    pub(crate) const TI23: f64 = 0.47662355450055045196;
    pub(crate) const TI31: f64 = -0.50287263494578687595;
    pub(crate) const TI32: f64 = 2.5719269498556054292;
    pub(crate) const TI33: f64 = -0.59603920482822492497;
    pub(crate) const C1: f64 = 0.1550510257216822;
    pub(crate) const C2: f64 = 0.6449489742783178;
    pub(crate) const C1M1: f64 = -0.8449489742783178;
    pub(crate) const C2M1: f64 = -0.3550510257216822;
    pub(crate) const C1MC2: f64 = -0.4898979485566356;
    pub(crate) const ALPHA: f64 = 2.681082873627752;
    pub(crate) const BETA: f64 = 3.050430199247411;
    pub(crate) const U1: f64 = 3.637834252744496;
}

impl Radau5 {
    pub fn integrate(
        &mut self,
        dydx: &'_ dyn Fn(ArrayViewMut1<f64>, ArrayView1<f64>, f64),
        jac: &'_ dyn Fn(ArrayViewMut2<f64>, ArrayView1<f64>, f64),
    ) -> (Vec<f64>, Vec<Array1<f64>>) {
        dbg!(&self.y);
        let mut xs: Vec<f64> = vec![];
        let mut ys: Vec<Array1<f64>> = vec![];

        let posneg = 1f64.copysign(self.xend - self.x);
        let hmaxn = self.hmax.abs().min((self.xend - self.x).abs());
        let cfac = self.safe * (1 + 2 * self.nit) as f64;

        self.h = self.h.abs().min(hmaxn).copysign(posneg);
        self.hold = self.h;

        let mut last = false;

        if (self.x + self.h * 1.0001 - self.xend) * posneg >= 0.0 {
            self.h = self.xend - self.x;
            last = true;
        }

        let mut hopt = self.h;
        let mut faccon = 1f64;

        for i in 0..self.n {
            self.cont[i] = self.y[i];
        }
        self.solution_output(&mut xs, &mut ys);

        for i in 0..self.n {
            self.scal[i] = self.atoler[i] + self.rtoler[i] * self.y[i].abs();
        }

        dydx(self.y0.view_mut(), self.y.view(), self.x);
        self.nfcn += 1;

        let mut hacc: f64 = 0.0;
        let mut erracc: f64 = 0.0;
        let mut thqold: f64 = 0.0;
        let mut nsing = 0;
        let mut ier;

        // basic integration step
        jac(self.fjac.view_mut(), self.y.view(), self.x);
        self.caljac = true;
        let mut loopp = true;
        while loopp {
            loopp = false;
            // compute the matrices e1 and e2 and their decompositions
            self.fac1 = Radau5::U1 / self.h;
            self.alphn = Radau5::ALPHA / self.h;
            self.betan = Radau5::BETA / self.h;

            ier = self.decomp_real();

            if ier != 0 {
                nsing += 1;
                if nsing >= 5 {
                    println!("Matrix is repeatedly singular");
                    return (xs, ys);
                }
                self.h *= 0.5;
                self.reject = true;
                last = false;
                if !self.caljac {
                    jac(self.fjac.view_mut(), self.y.view(), self.x);
                    self.caljac = true;
                }
                loopp = true;
                continue;
            }

            ier = self.decomp_complex();

            if ier != 0 {
                nsing += 1;
                if nsing >= 5 {
                    println!("Matrix is repeatedly singular");
                    return (xs, ys);
                }
                self.h *= 0.5;
                self.reject = true;
                last = false;
                if !self.caljac {
                    jac(self.fjac.view_mut(), self.y.view(), self.x);
                    self.caljac = true;
                }
                loopp = true;
                continue;
            }
            self.ndec += 1;

            loop {
                self.nstep += 1;
                if self.nstep >= self.nmax {
                    println!("More than {} steps needed.", self.nmax);
                    return (xs, ys);
                }

                if 0.1 * self.h.abs() <= self.x.abs() * f64::EPSILON {
                    println!("dx = {} less than dtmin.", self.h.abs());
                    return (xs, ys);
                }

                let xph = self.x + self.h;
                //  starting values for Newton iteration
                if self.first || self.startn {
                    for i in 0..self.n {
                        self.z1[i] = 0.0;
                        self.z2[i] = 0.0;
                        self.z3[i] = 0.0;
                        self.f1[i] = 0.0;
                        self.f2[i] = 0.0;
                        self.f3[i] = 0.0;
                    }
                } else {
                    let c3q = self.h / self.hold;
                    let c1q = Radau5::C1 * c3q;
                    let c2q = Radau5::C2 * c3q;
                    for i in 0..self.n {
                        let ak1 = self.cont[i + self.n];
                        let ak2 = self.cont[i + 2 * self.n];
                        let ak3 = self.cont[i + 3 * self.n];
                        self.z1[i] =
                            c1q * (ak1 + (c1q - Radau5::C2M1) * (ak2 + (c1q - Radau5::C1M1) * ak3));
                        self.z2[i] =
                            c2q * (ak1 + (c2q - Radau5::C2M1) * (ak2 + (c2q - Radau5::C1M1) * ak3));
                        self.z3[i] =
                            c3q * (ak1 + (c3q - Radau5::C2M1) * (ak2 + (c3q - Radau5::C1M1) * ak3));
                        self.f1[i] = Radau5::TI11 * self.z1[i]
                            + Radau5::TI12 * self.z2[i]
                            + Radau5::TI13 * self.z3[i];
                        self.f2[i] = Radau5::TI21 * self.z1[i]
                            + Radau5::TI22 * self.z2[i]
                            + Radau5::TI23 * self.z3[i];
                        self.f3[i] = Radau5::TI31 * self.z1[i]
                            + Radau5::TI32 * self.z2[i]
                            + Radau5::TI33 * self.z3[i];
                    }
                }

                //  loop for the simplified Newton iteration
                let mut newt = 0;
                faccon = faccon.max(f64::EPSILON).powf(0.8);
                let mut theta = self.thet.abs();
                let mut dyno;
                let mut dynold = 0.0;

                loop {
                    if newt >= self.nit {
                        if ier != 0 {
                            nsing += 1;
                            if nsing >= 5 {
                                println!("Matrix is repeatedly singular");
                                return (xs, ys);
                            }
                        }
                        self.h *= 0.5;
                        self.reject = true;
                        last = false;
                        if !self.caljac {
                            jac(self.fjac.view_mut(), self.y.view(), self.x);
                            self.caljac = true;
                        }
                        loopp = true;
                        break;
                    }
                    // compute the right-hand side
                    for i in 0..self.n {
                        self.cont[i] = self.y[i] + self.z1[i];
                    }
                    dydx(
                        self.z1.view_mut(),
                        self.cont.view(),
                        self.x + Radau5::C1 * self.h,
                    );

                    for i in 0..self.n {
                        self.cont[i] = self.y[i] + self.z2[i];
                    }
                    dydx(
                        self.z2.view_mut(),
                        self.cont.view(),
                        self.x + Radau5::C2 * self.h,
                    );

                    for i in 0..self.n {
                        self.cont[i] = self.y[i] + self.z3[i];
                    }

                    dydx(self.z3.view_mut(), self.cont.view(), xph);
                    self.nfcn += 3;

                    // solve the linear systems
                    for i in 0..self.n {
                        let a1 = self.z1[i];
                        let a2 = self.z2[i];
                        let a3 = self.z3[i];
                        self.z1[i] = Radau5::TI11 * a1 + Radau5::TI12 * a2 + Radau5::TI13 * a3;
                        self.z2[i] = Radau5::TI21 * a1 + Radau5::TI22 * a2 + Radau5::TI23 * a3;
                        self.z3[i] = Radau5::TI31 * a1 + Radau5::TI32 * a2 + Radau5::TI33 * a3;
                    }
                    self.linear_solve();
                    self.nsol += 1;
                    newt += 1;
                    dyno = 0.0;
                    let mut denom;
                    for i in 0..self.n {
                        denom = self.scal[i];
                        dyno = dyno
                            + (self.z1[i] / denom).powi(2)
                            + (self.z2[i] / denom).powi(2)
                            + (self.z3[i] / denom).powi(2);
                    }
                    dyno = (dyno / (3 * self.n) as f64).sqrt();
                    // bad convergence or number of iterations to large
                    if newt > 1 && newt < self.nit {
                        let thq = dyno / dynold;
                        if newt == 2 {
                            theta = thq;
                        } else {
                            theta = (thq * thqold).sqrt();
                        }
                        thqold = thq;
                        if theta < 0.99 {
                            faccon = theta / (1.0 - theta);
                            let dyth = faccon * dyno * theta.powi((self.nit - 1 - newt) as i32)
                                / self.fnewt;
                            if dyth >= 1.0 {
                                let qnewt = 1.0e-4f64.max(20.0f64.min(dyth));
                                let hhfac =
                                    0.8 * qnewt.powf(-1.0 / (4 + self.nit - 1 - newt) as f64);
                                self.h *= hhfac;
                                self.reject = true;
                                last = false;
                                if !self.caljac {
                                    jac(self.fjac.view_mut(), self.y.view(), self.x);
                                    self.caljac = true;
                                }
                                loopp = true;
                                break;
                            }
                        } else {
                            if ier != 0 {
                                nsing += 1;
                                if nsing >= 5 {
                                    println!("Matrix is repeatedly singular");
                                    return (xs, ys);
                                }
                            }
                            self.h *= 0.5;
                            self.reject = true;
                            last = false;
                            if !self.caljac {
                                jac(self.fjac.view_mut(), self.y.view(), self.x);
                                self.caljac = true;
                            }
                            loopp = true;
                            break;
                        }
                    }
                    dynold = dyno.max(f64::EPSILON);
                    for i in 0..self.n {
                        self.f1[i] = self.f1[i] + self.z1[i];
                        self.f2[i] = self.f2[i] + self.z2[i];
                        self.f3[i] = self.f3[i] + self.z3[i];
                        self.z1[i] = Radau5::T11 * self.f1[i]
                            + Radau5::T12 * self.f2[i]
                            + Radau5::T13 * self.f3[i];
                        self.z2[i] = Radau5::T21 * self.f1[i]
                            + Radau5::T22 * self.f2[i]
                            + Radau5::T23 * self.f3[i];
                        self.z3[i] = Radau5::T31 * self.f1[i] + self.f2[i];
                    }
                    if faccon * dyno <= self.fnewt {
                        break;
                    }
                }

                if loopp {
                    break;
                }

                // error estimation
                self.err = 0.0;
                self.error_estimate(&dydx);

                // computation of hnew -- require 0.2 <= hnew/h <= 8.
                let fac = self.safe.min(cfac / (newt + 2 * self.nit) as f64);
                let mut quot = self.facr.max(self.facl.min(self.err.powf(0.25) / fac));
                let mut hnew = self.h / quot;

                //  is the error small enough ?
                if self.err < 1.0 {
                    // step is accepted
                    self.first = false;
                    self.naccpt += 1;
                    if self.pred {
                        // predictive controller of Gustafsson
                        if self.naccpt > 1 {
                            let mut facgus = (hacc / self.h)
                                * (self.err * self.err / erracc).powf(0.25)
                                / self.safe;
                            facgus = self.facr.max(self.facl.min(facgus));
                            quot = quot.max(facgus);
                            hnew = self.h / quot;
                        }
                        hacc = self.h;
                        erracc = 1.0e-2f64.max(self.err);
                    }
                    self.xold = self.x;
                    self.hold = self.h;
                    self.x = xph;
                    for i in 0..self.n {
                        self.y[i] += self.z3[i];
                        self.cont[i + self.n] = (self.z2[i] - self.z3[i]) / Radau5::C2M1;
                        let ak = (self.z1[i] - self.z2[i]) / Radau5::C1MC2;
                        let mut acont3 = self.z1[i] / Radau5::C1;
                        acont3 = (ak - acont3) / Radau5::C2;
                        self.cont[i + 2 * self.n] = (ak - self.cont[i + self.n]) / Radau5::C1M1;
                        self.cont[i + 3 * self.n] = self.cont[i + 2 * self.n] - acont3;
                    }
                    for i in 0..self.n {
                        self.scal[i] = self.atoler[i] + self.rtoler[i] * self.y[i].abs();
                    }
                    self.solution_output(&mut xs, &mut ys);
                    self.caljac = false;
                    if last {
                        self.h = hopt;
                        return (xs, ys);
                    }
                    dydx(self.y0.view_mut(), self.y.view(), self.x);
                    self.nfcn += 1;
                    hnew = posneg * hnew.abs().min(hmaxn);
                    hopt = self.h.min(hnew);
                    if self.reject {
                        hnew = posneg * hnew.abs().min(self.h.abs());
                    }
                    self.reject = false;
                    if (self.x + hnew / self.quot1 - self.xend) * posneg >= 0.0 {
                        self.h = self.xend - self.x;
                        last = true;
                    } else {
                        let qt = hnew / self.h;
                        if (theta <= self.thet) && (qt >= self.quot1) && (qt <= self.quot2) {
                            continue;
                        }
                        self.h = hnew;
                    }
                    if theta > self.thet {
                        jac(self.fjac.view_mut(), self.y.view(), self.x);
                        self.caljac = true;
                    }
                    loopp = true;
                } else {
                    // step is rejected
                    self.reject = true;
                    last = false;
                    if self.first {
                        self.h *= 0.1;
                    } else {
                        self.h = hnew;
                    }
                    if self.naccpt >= 1 {
                        self.nrejct += 1;
                    }
                    if !self.caljac {
                        jac(self.fjac.view_mut(), self.y.view(), self.x);
                        self.caljac = true;
                    }
                    loopp = true;
                }
                break;
            }
        }
        unreachable!()
    }
    fn decomp_real(&mut self) -> usize {
        // mass = identity, Jacobian a full matrix
        for j in 0..self.n {
            for i in 0..self.n {
                self.e1[[i, j]] = -self.fjac[[i, j]];
            }
            self.e1[[j, j]] += self.fac1;
        }
        dec(self.n, self.e1.view_mut(), self.ip1.view_mut())
    }
    fn decomp_complex(&mut self) -> usize {
        for j in 0..self.n {
            for i in 0..self.n {
                self.e2r[[i, j]] = -self.fjac[[i, j]];
                self.e2i[[i, j]] = 0.0;
            }
            self.e2r[[j, j]] += self.alphn;
            self.e2i[[j, j]] = self.betan;
        }
        decc(
            self.n,
            self.e2r.view_mut(),
            self.e2i.view_mut(),
            self.ip2.view_mut(),
        )
    }
    fn linear_solve(&mut self) {
        for i in 0..self.n {
            let s2 = -self.f2[i];
            let s3 = -self.f3[i];
            self.z1[i] -= self.f1[i] * self.fac1;
            self.z2[i] = self.z2[i] + s2 * self.alphn - s3 * self.betan;
            self.z3[i] = self.z3[i] + s3 * self.alphn + s2 * self.betan;
        }
        sol(self.n, self.e1.view(), self.z1.view_mut(), self.ip1.view());
        solc(
            self.n,
            self.e2r.view(),
            self.e2i.view(),
            self.z2.view_mut(),
            self.z3.view_mut(),
            self.ip2.view(),
        );
    }
    fn error_estimate(&mut self, dydx: &'_ &dyn Fn(ArrayViewMut1<f64>, ArrayView1<f64>, f64)) {
        let hee1 = -(13.0 + 7.0 * 6.0f64.sqrt()) / (3.0 * self.h);
        let hee2 = (-13.0 + 7.0 * 6.0f64.sqrt()) / (3.0 * self.h);
        let hee3 = -1.0 / (3.0 * self.h);

        for i in 0..self.n {
            self.f2[i] = hee1 * self.z1[i] + hee2 * self.z2[i] + hee3 * self.z3[i];
            self.cont[i] = self.f2[i] + self.y0[i];
        }
        sol(
            self.n,
            self.e1.view(),
            self.cont.view_mut(),
            self.ip1.view(),
        );

        self.err = 0.0;
        for i in 0..self.n {
            self.err += (self.cont[i] / self.scal[i]).powi(2);
        }
        self.err = (self.err / self.n as f64).sqrt().max(1.0e-10);

        if self.err < 1.0 {
            return;
        }

        if self.first || self.reject {
            for i in 0..self.n {
                self.cont[i] += self.y[i];
            }
            dydx(self.f1.view_mut(), self.cont.view(), self.x);
            for i in 0..self.n {
                self.cont[i] = self.f1[i] + self.f2[i];
            }

            sol(
                self.n,
                self.e1.view(),
                self.cont.view_mut(),
                self.ip1.view(),
            );

            self.err = 0.0;
            for i in 0..self.n {
                self.err += (self.cont[i] / self.scal[i]).powi(2);
            }
            self.err = (self.err / self.n as f64).sqrt().max(1.0e-10);
        }
    }
    fn continuous_output(&self) -> Array1<f64> {
        let sq6 = 6f64.sqrt();
        let c1 = (4.0 - sq6) / 10.0;
        let c2 = (4.0 + sq6) / 10.0;
        let c1m1 = c1 - 1.0;
        let c2m1 = c2 - 1.0;

        let s = (self.xd - self.x) / self.hold;
        let mut out = Array1::zeros(self.n);
        for i in 0..self.n {
            out[i] = self.cont[i]
                + s * (self.cont[i + self.n]
                    + (s - c2m1)
                        * (self.cont[i + 2 * self.n] + (s - c1m1) * self.cont[i + 3 * self.n]));
        }
        out
    }
    fn solution_output(&mut self, xs: &mut Vec<f64>, ys: &mut Vec<Array1<f64>>) {
        if self.naccpt == 0 {
            self.xd = self.xold;
        }

        while self.xd < self.x {
            if self.xold <= self.xd && self.x >= self.xd {
                xs.push(self.xd);
                ys.push(self.continuous_output());
                self.xd += self.dx;
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_van_der_pol() {
        let dydx = |mut dy: ArrayViewMut1<f64>, y: ArrayView1<f64>, x: f64| {
            dy[0] = y[1];
            dy[1] = (1.0 - y[0] * y[0]) * y[1] - y[0];
        };
        let jac = |mut j: ArrayViewMut2<f64>, y: ArrayView1<f64>, x: f64| {
            j[[0, 0]] = 0.0;
            j[[0, 1]] = 1.0;
            j[[1, 0]] = -2.0 * y[0] * y[1] - 1.0;
            j[[1, 1]] = 1.0 - y[0] * y[0];
        };
        let y: Array1<f64> = Array::from(vec![1.0, 0.0]);
        let mut integrator = Radau5::builder(y, (0.0, 10.0)).dx(1e-2).build().unwrap();
        let (xs, ys) = integrator.integrate(&dydx, &jac);

        for (x, y) in xs.iter().zip(ys.iter()) {
            println!("x, y = {}, {:?}", x, y);
        }
    }
}
