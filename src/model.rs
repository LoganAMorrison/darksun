//pub mod boltz;
pub mod cross_sections;
pub mod data;
pub mod thermal_functions;

#[cfg(test)]
mod test {
    use super::*;
    use crate::model::thermal_functions::*;
    use crate::ode::radau5::*;
    use crate::standard_model::*;
    use crate::DarkSun;
    use haliax_constants::cosmology::T_CMB;
    use ndarray::prelude::*;

    #[test]
    fn test_boltzmann() {
        let mut model = DarkSun::new(6, 1e-3);
        let td = model.lamc / 2.0;
        let xi = model.compute_xi_const_td(td).unwrap();
        let tsm = td / xi;

        let xstart = model.m_eta / tsm;
        let xfinal = model.m_eta / T_CMB;
        let s = sm_entropy_density(tsm);
        let w_eta = (model.neq_eta(td) / s).ln();
        let y_del = (model.neq_del(td) / s) * (-model.adel * model.n as f64).exp();

        let logx_span = (xstart.ln(), xfinal.ln());
        let dlogx = (logx_span.1 - logx_span.0) / 100.0;

        let w = Array::from(vec![w_eta, y_del]);
        let mut rad = Radau5::builder(w, logx_span).dx(dlogx).build().unwrap();
        let (logxs, ws) = rad.integrate(&mut model).unwrap();
        for (logx, w) in logxs.iter().zip(ws.iter()) {
            println!("{}, {:e}, {:e}", logx, w[0].exp(), w[1]);
        }
        println!("xifo = {:?}", model.xi_fo);
        println!(
            "tbfo = {:?}",
            model.m_eta / (model.xi_fo.unwrap() * model.tsm_fo.unwrap())
        );
        println!("xibb = {:?}", model.xi_bbn);
    }
}
