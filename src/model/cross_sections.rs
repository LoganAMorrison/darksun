use super::data::*;
use crate::DarkSun;
use cyphus_integration::prelude::*;
use cyphus_specfun::bessel::CylBesselK;
use std::f64::consts::PI;

impl DarkSun {
    /// Cross section for 2eta -> 4eta with all prefactors included as a
    /// function of the center-of-mass energy.
    pub fn rust_cross_section_2eta_4eta(&self, cme: f64) -> f64 {
        let mu = self.mu_eta;
        let n = self.n as f64;
        let lam = self.lamc;
        //         v = (256 pi^4/9)^2
        let norm = 6.9093374296577904e7 * mu.powi(14) / (lam * lam * n.powi(11));
        let c44 = norm * self.lec1.powi(4) / 9.0;
        let c66 = norm * self.lec2.powi(2) / 25.0;
        let c46 = -2.0 * norm * self.lec2 * self.lec1.powi(2) / 15.0;

        let z = cme / self.m_eta;
        c44 * scaled_cs_eta_24_44(z) + c46 * scaled_cs_eta_24_46(z) + c66 * scaled_cs_eta_24_66(z)
    }
    /// Compute the thermally-averaged annihilation cross section for
    /// 2eta -> 4eta as a function of x = meta/Teta,
    pub fn rust_thermal_cross_section_2eta_4eta(&self, x: f64) -> f64 {
        let den = 2.0 * x.cyl_bessel_kn_scaled(2);
        let pre = x / (den * den);

        let f = |z: f64| {
            let z2 = z * z;
            let sig = self.rust_cross_section_2eta_4eta(z * self.m_eta);
            let kernal = z2 * (z2 - 4.0) * (x * z).cyl_bessel_k1_scaled() * (-x * (z - 2.0)).exp();
            sig * kernal
        };
        let int = qagi(f, 4.0, f64::INFINITY, 1e-8, 1e-8, 500, 2).val;
        int * pre
    }
    /// Compute the thermally-averaged annihilation cross section for
    /// 4eta -> 2eta as a function of x = meta/Teta,
    pub fn rust_thermal_cross_section_4eta_2eta(&self, x: f64) -> f64 {
        let bes = x.cyl_bessel_kn_scaled(2);
        let pre = PI.powi(4) * x.powi(3) / (self.m_eta.powi(6) * bes.powi(4));

        let f = |z: f64| {
            let z2 = z * z;
            let sig = self.rust_cross_section_2eta_4eta(z * self.m_eta);
            let kernal = z2 * (z2 - 4.0) * (x * z).cyl_bessel_k1_scaled() * (-x * (z - 4.0)).exp();
            sig * kernal
        };
        let int = qagi(f, 4.0, f64::INFINITY, 1e-8, 1e-8, 500, 2).val;
        int * pre
    }
    /// Compute the thermally-averaged annihilation cross section for
    /// 2eta -> 2delta as a function of x = meta/Teta,
    pub fn rust_thermal_cross_section_2eta_2del(&self, x: f64) -> f64 {
        let c = self.c;
        let n = self.n as f64;
        let lam = self.lamc;
        let den = 2.0 * x.cyl_bessel_kn_scaled(2);
        let pf = x / (den * den);
        let zmin = 2.0 * self.m_del / self.m_eta;
        let sig = (-2.0 * c * n).exp() / (64.0 * PI * n * n * lam * lam);

        let f = |z: f64| {
            let z2 = z * z;
            let kernal = z2 * (z2 - 4.0) * (x * z).cyl_bessel_k1_scaled() * (-x * (z - 2.0)).exp();
            kernal
        };
        let int = qagi(f, zmin, f64::INFINITY, 1e-8, 1e-8, 500, 2).val;

        pf * sig * int
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_thermal_cs() {
        let model = DarkSun::new(20, 1e-3);
        model.thermal_cross_section_2eta_4eta(1.0);
    }
}
