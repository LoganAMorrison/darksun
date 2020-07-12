//! The Dark SU(N) Model.
//! This file contains the python API to the dark SU(N) model.
pub mod brent;
pub mod interp_utils;
pub mod model;
pub mod ode;
pub mod standard_model;

use crate::ode::radau5::*;
use crate::standard_model::*;
use haliax_constants::cosmology::*;
use ndarray::prelude::*;
use pyo3::exceptions;
use pyo3::prelude::*;

#[pymodule]
/// The Dark Sun module. This module contains the `DarkSun` class which has
/// all the functionalitiy to study the Dark SU(N) model.
fn darksun(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_class::<DarkSun>()?;
    Ok(())
}

#[pyclass(name = DarkSun, module = "darksun", dict)]
/// Class for computing all aspects of the Dark SU(N) model.
pub struct DarkSun {
    /// Order of the SU(N) gauge group.
    pub n: usize,
    /// Confinement scale.
    pub lamc: f64,
    /// Temperature ration td/tsm above confinement scale.
    #[pyo3(get, set)]
    pub xi_inf: f64,
    /// Constant parameterizing delta+delta->eta+eta matrix element. It is
    /// defined such that M ~ exp(-c*n)
    #[pyo3(get, set)]
    pub c: f64,
    /// Low-energy constants in 4-pt eta interaction
    #[pyo3(get, set)]
    pub lec1: f64,
    /// Low-energy constant in the 6-pt eta interaction
    #[pyo3(get, set)]
    pub lec2: f64,
    /// Exponential supression factor for initial delta density
    #[pyo3(get, set)]
    pub adel: f64,
    /// Prefactor of the eta mass defined such that
    /// meta = mu_eta * lamc / sqrt(n).
    pub mu_eta: f64,
    /// Prefactor of the delta mass defined such that
    /// mdel = mu_del * lamc * n.
    pub mu_del: f64,

    /// Mass of the eta.
    #[pyo3(get)]
    pub m_eta: f64,
    /// Mass of the delta.
    #[pyo3(get)]
    pub m_del: f64,
    #[pyo3(get)]
    /// Value of xi = Td / Tsm when the eta freezes out.
    pub xi_fo: Option<f64>,
    #[pyo3(get)]
    /// Value of Tsm when the eta freezes out.
    pub tsm_fo: Option<f64>,
    #[pyo3(get)]
    /// Value of xi = Td / Tsm at BBN.
    pub xi_bbn: Option<f64>,
    #[pyo3(get)]
    /// Value of xi = Td / Tsm at CMB.
    pub xi_cmb: Option<f64>,
    #[pyo3(get)]
    /// Relic density of the eta.
    pub rd_eta: Option<f64>,
    #[pyo3(get)]
    /// Relic density of the delta.
    pub rd_del: Option<f64>,
    #[pyo3(get)]
    /// dNeff at CMB.
    pub dneff_cmb: Option<f64>,
    #[pyo3(get)]
    /// dNeff at BBN.
    pub dneff_bbn: Option<f64>,
    #[pyo3(get)]
    /// Value of sigma_eta / m_eta.
    pub eta_si_per_mass: Option<f64>,
    #[pyo3(get)]
    /// Value of sigma_del / m_del.
    pub del_si_per_mass: Option<f64>,
}

#[pymethods]
impl DarkSun {
    #[new]
    /// Create a new DarkSun object
    pub fn new(n: usize, lamc: f64) -> Self {
        DarkSun {
            n: n,
            lamc: lamc,
            xi_inf: 1e-2,
            c: 1.0,
            lec1: 1.0,
            lec2: 1.0,
            adel: 1.0,
            mu_eta: 1.0,
            mu_del: 1.0,
            m_eta: lamc / (n as f64).sqrt(),
            m_del: lamc * n as f64,
            xi_fo: None,
            tsm_fo: None,
            xi_bbn: None,
            xi_cmb: None,
            rd_eta: None,
            rd_del: None,
            dneff_cmb: None,
            dneff_bbn: None,
            eta_si_per_mass: None,
            del_si_per_mass: None,
        }
    }
    #[getter]
    pub fn get_n(&self) -> PyResult<usize> {
        Ok(self.n)
    }
    #[getter]
    pub fn get_lamc(&self) -> PyResult<f64> {
        Ok(self.lamc)
    }
    #[getter]
    pub fn get_mu_eta(&self) -> PyResult<f64> {
        Ok(self.mu_eta)
    }
    #[getter]
    pub fn get_mu_del(&self) -> PyResult<f64> {
        Ok(self.mu_del)
    }

    #[setter]
    pub fn set_n(&mut self, n: usize) -> PyResult<()> {
        self.n = n;
        self.m_eta = self.mu_eta * self.lamc / (self.n as f64).sqrt();
        self.m_del = self.mu_del * self.lamc * (self.n as f64);
        Ok(())
    }
    #[setter]
    pub fn set_lamc(&mut self, lamc: f64) -> PyResult<()> {
        self.lamc = lamc;
        self.m_eta = self.mu_eta * self.lamc / (self.n as f64).sqrt();
        self.m_del = self.mu_del * self.lamc * (self.n as f64);
        Ok(())
    }
    #[setter]
    pub fn set_mu_eta(&mut self, mu_eta: f64) -> PyResult<()> {
        self.mu_eta = mu_eta;
        self.m_eta = self.mu_eta * self.lamc / (self.n as f64).sqrt();
        self.m_del = self.mu_del * self.lamc * (self.n as f64);
        Ok(())
    }
    #[setter]
    pub fn set_mu_del(&mut self, mu_del: f64) -> PyResult<()> {
        self.mu_del = mu_del;
        self.m_eta = self.mu_eta * self.lamc / (self.n as f64).sqrt();
        self.m_del = self.mu_del * self.lamc * (self.n as f64);
        Ok(())
    }
    /// Compute the value of xi = Td / Tsm given the standard model temperature.
    /// This computation is done by using a root finding algorithm to solve
    /// a non-linear equation. A `ValueError` is raised if the root-finding
    /// algorithm fails.
    pub fn compute_xi(&self, tsm: f64) -> PyResult<f64> {
        let res = match self.tsm_fo {
            Some(tfo) => {
                let xi_fo = self.xi_fo.unwrap();
                if tfo * xi_fo > self.m_eta {
                    Ok(xi_fo)
                } else {
                    Ok(xi_fo * tsm / tfo)
                }
            }
            None => self.compute_xi_const_tsm(tsm),
        };
        match res {
            Ok(val) => Ok(val),
            Err(err) => Err(exceptions::ValueError::py_err(format!(
                "Failed to compute xi: {:?}",
                err
            ))),
        }
    }
    /// Cross section for 2eta -> 4eta with all prefactors included as a
    /// function of the center-of-mass energy.
    pub fn cross_section_2eta_4eta(&self, cme: f64) -> PyResult<f64> {
        Ok(self.rust_cross_section_2eta_4eta(cme))
    }
    /// Compute the thermally-averaged annihilation cross section for
    /// 2eta -> 4eta as a function of x = meta/Teta,
    pub fn thermal_cross_section_2eta_4eta(&self, x: f64) -> PyResult<f64> {
        Ok(self.rust_thermal_cross_section_2eta_4eta(x))
    }
    /// Compute the thermally-averaged annihilation cross section for
    /// 2eta -> 4eta as a function of x = meta/Teta,
    pub fn thermal_cross_section_4eta_2eta(&self, x: f64) -> PyResult<f64> {
        Ok(self.rust_thermal_cross_section_4eta_2eta(x))
    }
    /// Compute the thermally-averaged annihilation cross section for
    /// 2eta -> 2delta as a function of x = meta/Teta,
    pub fn thermal_cross_section_2eta_2delta(&self, x: f64) -> PyResult<f64> {
        Ok(self.rust_thermal_cross_section_2eta_2del(x))
    }
    pub fn solve_boltzmann(&mut self) -> PyResult<(Vec<f64>, Vec<Vec<f64>>)> {
        let td = self.lamc / 2.0;
        let xi = self.compute_xi_const_td(td).unwrap();
        let tsm = td / xi;

        let xstart = self.m_eta / tsm;
        let xfinal = self.m_eta / T_CMB;
        let s = sm_entropy_density(tsm);
        let w_eta = (self.neq_eta(td) / s).ln();
        let y_del = (self.neq_del(td) / s) * (-self.adel * self.n as f64).exp();

        let logx_span = (xstart.ln(), xfinal.ln());
        let dlogx = (logx_span.1 - logx_span.0) / 100.0;

        let w = Array::from(vec![w_eta, y_del]);
        let mut rad = Radau5::builder(w, logx_span).dx(dlogx).build();
        match rad {
            Ok(mut integrator) => {
                let solution = integrator.integrate(self);
                match solution {
                    Ok(sol) => {
                        let mut us: Vec<Vec<f64>> = vec![];
                        us.resize(sol.1.len(), vec![0.0; (sol.1)[0].len()]);
                        for i in 0..us.len() {
                            us[i] = (sol.1)[i].to_vec();
                        }
                        Ok((sol.0, us))
                    }
                    Err(err) => Err(exceptions::RuntimeError::py_err(err)),
                }
            }
            Err(err) => Err(exceptions::RuntimeError::py_err(err)),
        }
    }
}
