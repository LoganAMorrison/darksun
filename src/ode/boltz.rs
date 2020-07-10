use crate::DarkSun;
use haliax_constants::prelude::*;
use haliax_thermal_functions::prelude::*;
use ndarray::prelude::*;
use std::f64::consts::PI;

impl DarkSun {
    fn sqrt_gstar(&self, tsm: f64, xi: f64) -> f64 {
        let td = xi * tsm;
        let gd = self.dark_geff(td);
        let gsm = sm_geff(tsm);
        return sm_sqrt_gstar(tsm) * (gsm / (gsm + gd * xi * xi * xi * xi)).sqrt();
    }
    pub fn dudt(&self, mut dw: ArrayViewMut1<f64>, w: ArrayView1<f64>, logx: f64) {
        let x = logx.exp();

        let meta = self.m_eta;
        let tsm = meta / x;
        let xi = self.compute_xi(tsm).unwrap();
        let td = xi * tsm;
        let s = sm_entropy_density(tsm);
        let neq = self.neq_eta(td);
        let weq = (neq / s).ln();

        let sig_ee_eeee = self.rust_thermal_cross_section_2eta_4eta(meta / td);
        let sig_eeee_ee = sig_ee_eeee / neq / neq;

        let pf_e = -(PI / 45.0).sqrt() * M_PLANK * self.sqrt_gstar(tsm, xi) * tsm * s * s;

        dw[0] = pf_e * sig_eeee_ee * w[0].exp() * ((2.0 * w[0]).exp() - (2.0 * weq).exp());
    }

    pub fn dfdu(&mut self, mut j: ArrayViewMut2<f64>, w: ArrayView1<f64>, logx: f64) {
        let x = logx.exp();

        let meta = self.m_eta;
        let tsm = meta / x;
        let xi = self.compute_xi(tsm).unwrap();
        let td = xi * tsm;
        let s = sm_entropy_density(tsm);
        let neq = self.neq_eta(td);
        let weq = (neq / s).ln();

        let sig_ee_eeee = self.rust_thermal_cross_section_2eta_4eta(meta / td);
        let sig_eeee_ee = sig_ee_eeee / neq / neq;

        let pf_e = -(PI / 45.0).sqrt() * M_PLANK * self.sqrt_gstar(tsm, xi) * tsm * s * s;
        j[[0, 0]] =
            pf_e * sig_eeee_ee * w[0].exp() * (3.0 * (2.0 * w[0]).exp() - (2.0 * weq).exp());
    }
}
