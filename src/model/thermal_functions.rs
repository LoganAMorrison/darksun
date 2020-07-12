use crate::brent::brent;
use crate::standard_model::*;
use crate::DarkSun;
use cyphus_specfun::bessel::CylBesselK;
use cyphus_specfun::lambert_w::LambertW;
use std::f64::consts::PI;

const SM_HEFF_INF: f64 = 106.83;
const SM_HEFF_0: f64 = 3.93879;

impl DarkSun {
    pub fn neq_eta(&self, td: f64) -> f64 {
        let x = self.m_eta / td;
        td * td * td / (2.0 * PI * PI) * x * x * x.cyl_bessel_kn(2)
    }
    pub fn neq_del(&self, td: f64) -> f64 {
        let x = self.m_del / td;
        let g = (self.n + 1) as f64;
        g * td * td * td / (2.0 * PI * PI) * x * x * x.cyl_bessel_kn(2)
    }
    pub fn heff_eta(&self, td: f64) -> f64 {
        let x = self.m_eta / td;
        let x3 = x * x * x;
        let pi4 = PI * PI * PI * PI;
        45.0 / (4.0 * pi4) * x3 * x.cyl_bessel_kn(3)
    }
    pub fn heff_del(&self, td: f64) -> f64 {
        let x = self.m_del / td;
        let x3 = x * x * x;
        let pi4 = PI * PI * PI * PI;
        let g = (self.n + 1) as f64;
        g * 45.0 / (4.0 * pi4) * x3 * x.cyl_bessel_kn(3)
    }
    pub fn dark_heff(&self, td: f64) -> f64 {
        let pi4 = PI * PI * PI * PI;
        let xe = self.m_eta / td;
        let xe3 = xe * xe * xe;
        let xd = self.m_del / td;
        let gd = (self.n + 1) as f64;
        let xd3 = xd * xd * xd;
        45.0 / (4.0 * pi4) * (xe3 * xe.cyl_bessel_kn(3) + gd * xd3 * xd.cyl_bessel_kn(3))
    }
    pub fn dark_geff(&self, td: f64) -> f64 {
        let pi4 = PI * PI * PI * PI;
        let xe = self.m_eta / td;
        let xd = self.m_del / td;
        let gd = (self.n + 1) as f64;
        30.0 / (2.0 * pi4)
            * (xe * xe * (xe * xe.cyl_bessel_k1() + 3.0 * xe.cyl_bessel_kn(2))
                + gd * xd * xd * (xd * xd.cyl_bessel_kn(3) + 3.0 * xd.cyl_bessel_kn(2)))
    }
}

impl DarkSun {
    fn hd_inf(&self) -> f64 {
        let n = self.n as f64;
        7.0 / 2.0 * n + 2.0 * n * n
    }
    fn sum_g(&self) -> f64 {
        1.0 + (self.n as f64 + 1.0)
    }
    fn gl(&self) -> f64 {
        1.0
    }
    fn ml(&self) -> f64 {
        self.m_eta
    }

    /// Compute the upper bound on xi assuming td is known.
    fn xi_upper_bound_const_td(&self, td: f64) -> f64 {
        let hd = self.dark_heff(td);
        (self.hd_inf() / hd * SM_HEFF_0 / SM_HEFF_INF).cbrt() * self.xi_inf
    }
    ///Compute the lower bound on xi assuming td is known.
    fn xi_lower_bound_const_td(&self, td: f64) -> f64 {
        let hd = self.dark_heff(td);
        (self.hd_inf() / hd).cbrt() * self.xi_inf
    }

    /// Compute the upper bound on xi assuming tsm is known.
    fn xi_lower_bound_const_tsm(&self, tsm: f64) -> f64 {
        (sm_heff(tsm) * self.hd_inf() / self.sum_g() / SM_HEFF_INF).cbrt() * self.xi_inf
    }

    /// Compute the lower bound on xi assuming tsm is known.
    fn xi_upper_bound_const_tsm(&self, tsm: f64) -> f64 {
        let xl = self.m_eta / tsm;
        let hsm = sm_heff(tsm);

        let lw_arg_num = (45.0 * SM_HEFF_INF * xl * xl * xl).powi(2);
        let lw_arg_den = (4.0 * self.hd_inf() * hsm * self.xi_inf).powi(2) * PI.powi(7);
        return 2.0 * xl / (lw_arg_num / lw_arg_den).lambert_w0();
    }
    /// Compute xi = Td/Tsm assuming td is known
    pub fn compute_xi_const_td(&self, td: f64) -> Result<f64, String> {
        let hd = self.dark_heff(td);
        let c1 = self.hd_inf() * self.xi_inf.powi(3) / SM_HEFF_INF;
        let lb = 0.8 * self.xi_lower_bound_const_td(td);
        let ub = 1.2 * self.xi_upper_bound_const_td(td);

        let f = |xi: f64| hd * xi.powi(3) - sm_heff(td / xi) * c1;

        brent(f, lb, ub, 1e-5)
    }
    /// Compute xi = Td/Tsm assuming tsm is known
    pub fn compute_xi_const_tsm(&self, tsm: f64) -> Result<f64, String> {
        match self.tsm_fo {
            Some(tfo) => {
                let xi_fo = self.xi_fo.unwrap();
                if tfo * xi_fo > self.m_eta {
                    Ok(xi_fo)
                } else {
                    Ok(xi_fo * tsm / tfo)
                }
            }
            None => {
                let ub = 1.2 * self.xi_upper_bound_const_tsm(tsm);
                let lb = 0.8 * self.xi_lower_bound_const_tsm(tsm);
                let hsm = sm_heff(tsm);
                let c1 = hsm * self.hd_inf() * self.xi_inf.powi(3) / SM_HEFF_INF;

                let f = |xi: f64| self.dark_heff(xi * tsm) * xi.powi(3) - c1;

                brent(f, lb, ub, 1e-5)
            }
        }
    }
}
