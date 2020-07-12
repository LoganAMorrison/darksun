use crate::standard_model::*;
use crate::DarkSun;
use haliax_constants::prelude::*;
use ndarray::prelude::*;
use std::f64::consts::PI;

impl DarkSun {
    fn sqrt_gstar(&self, tsm: f64, xi: f64) -> f64 {
        let td = xi * tsm;
        let gd = self.dark_geff(td);
        let gsm = sm_geff(tsm);
        return sm_sqrt_gstar(tsm) * (gsm / (gsm + gd * xi * xi * xi * xi)).sqrt();
    }
    pub fn boltzmann(&mut self, mut dw: ArrayViewMut1<f64>, w: ArrayView1<f64>, logx: f64) {
        let x = logx.exp();

        let meta = self.m_eta;
        let tsm = meta / x;
        let xi = self.compute_xi(tsm).unwrap();
        let td = xi * tsm;
        let s = sm_entropy_density(tsm);
        let neq = self.neq_eta(td);
        let weq = (neq / s).ln();

        if w[0] - weq > 0.1 && self.xi_fo.is_none() {
            self.xi_fo = Some(xi);
            self.tsm_fo = Some(tsm);
        }
        if tsm > T_BBN && self.xi_bbn.is_none() {
            self.xi_bbn = Some(xi);
        }

        let sige = self.rust_thermal_cross_section_4eta_2eta(meta / td);
        let sigd = self.rust_thermal_cross_section_2eta_2del(meta / td);

        let fac = (PI / 45.0).sqrt() * M_PLANK * self.sqrt_gstar(tsm, xi);
        let pfe = -fac * tsm * s * s;
        let pfd = fac * meta / (x * x);

        dw[0] = pfe * sige * w[0].exp() * ((2.0 * w[0]).exp() - (2.0 * weq).exp());
        dw[1] = pfd * sigd * (2.0 * w[0]).exp();
    }
    pub fn boltzmann_jac(&self, mut j: ArrayViewMut2<f64>, w: ArrayView1<f64>, logx: f64) {
        let x = logx.exp();

        let meta = self.m_eta;
        let tsm = meta / x;
        let xi = self.compute_xi(tsm).unwrap();
        let td = xi * tsm;
        let s = sm_entropy_density(tsm);
        let neq = self.neq_eta(td);
        let weq = (neq / s).ln();

        let sige = self.rust_thermal_cross_section_4eta_2eta(meta / td);
        let sigd = self.rust_thermal_cross_section_2eta_2del(meta / td);

        let fac = (PI / 45.0).sqrt() * M_PLANK * self.sqrt_gstar(tsm, xi);
        let pfe = -fac * tsm * s * s;
        let pfd = fac * meta / (x * x);

        j[[0, 0]] = pfe * sige * w[0].exp() * (3.0 * (2.0 * w[0]).exp() - (2.0 * weq).exp());
        j[[0, 1]] = 0.0;
        j[[1, 0]] = 2.0 * pfd * sigd * (2.0 * w[0]).exp();
        j[[1, 1]] = 0.0;
    }
}
