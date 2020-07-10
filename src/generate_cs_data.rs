//! File for generating the cross-section data for the integrated squared matrix
//! element for 2eta->4eta.

use haliax_phase_space::prelude::*;
use indicatif::ProgressBar;
use ndarray::prelude::*;
use std::io::prelude::*;

fn main() -> std::io::Result<()> {
    let mut file_zs = std::fs::File::create("data/log_scaled_cs_data_zs.dat")?;
    let mut file_cs44 = std::fs::File::create("data/log_scaled_cs_data_cs44.dat")?;
    let mut file_cs66 = std::fs::File::create("data/log_scaled_cs_data_cs66.dat")?;
    let mut file_cs46 = std::fs::File::create("data/log_scaled_cs_data_cs46.dat")?;

    let zs = Array1::<f64>::logspace(10.0, (4f64 + 1e-5).log10(), 2.0, 500);
    let bar = ProgressBar::new(zs.len() as u64);

    for z in zs.iter() {
        let cs = scaled_cross_section_2_4(*z);
        //        println!(
        //            "z, c44, c66, c46, tot = {:e}, {:e}, {:e}, {:e}",
        //            z,
        //            cs.0,
        //            cs.1,
        //            cs.2,
        //        );
        file_zs.write(format!("{}\n", z.log10()).to_string().as_bytes())?;
        file_cs44.write(format!("{}\n", cs.0.log10()).to_string().as_bytes())?;
        file_cs66.write(format!("{}\n", cs.1.log10()).to_string().as_bytes())?;
        file_cs46.write(format!("{}\n", cs.2.log10()).to_string().as_bytes())?;
        bar.inc(1);
    }
    bar.finish();
    Ok(())
}

/// Template matrix element squared for the various diagrams that contribute
/// to 2eta <-> 4eta using only 4pt interactions.
fn mat_elem_2_to_4_diag_4pt(
    p1: &FourMomentum,
    p2: &FourMomentum,
    p3: &FourMomentum,
    p4: &FourMomentum,
    p5: &FourMomentum,
    p6: &FourMomentum,
    q: &FourMomentum,
) -> f64 {
    (p1.dot(&q) * p2.dot(&p3) + p2.dot(&q) * p3.dot(&p1) + p3.dot(&q) * p1.dot(&p2))
        * (p4.dot(&q) * p5.dot(&p6) + p5.dot(&q) * p6.dot(&p4) + p6.dot(&q) * p4.dot(&p5))
        / (q.dot(&q) - 1.0)
}
/// Template matrix element squared for the various diagrams that contribute
/// to 2eta <-> 4eta using only 6pt interactions.
#[allow(dead_code)]
fn mat_elem_2_to_4_diag_6pt(
    p1: &FourMomentum,
    p2: &FourMomentum,
    p3: &FourMomentum,
    p4: &FourMomentum,
    p5: &FourMomentum,
    p6: &FourMomentum,
) -> f64 {
    p1.dot(&p4) * p2.dot(&p6) * p3.dot(&p5)
        + p1.dot(&p4) * p2.dot(&p5) * p3.dot(&p6)
        + p1.dot(&p3) * p2.dot(&p6) * p4.dot(&p5)
        + p1.dot(&p2) * p3.dot(&p6) * p4.dot(&p5)
        + p1.dot(&p6)
            * (p2.dot(&p5) * p3.dot(&p4) + p2.dot(&p4) * p3.dot(&p5) + p2.dot(&p3) * p4.dot(&p5))
        + p1.dot(&p3) * p2.dot(&p5) * p4.dot(&p6)
        + p1.dot(&p2) * p3.dot(&p5) * p4.dot(&p6)
        + p1.dot(&p5)
            * (p2.dot(&p6) * p3.dot(&p4) + p2.dot(&p4) * p3.dot(&p6) + p2.dot(&p3) * p4.dot(&p6))
        + p1.dot(&p4) * p2.dot(&p3) * p5.dot(&p6)
        + p1.dot(&p3) * p2.dot(&p4) * p5.dot(&p6)
        + p1.dot(&p2) * p3.dot(&p4) * p5.dot(&p6)
}
/// Compute the integrated squared matrix element for 2eta -> 4eta in the
/// center of mass frame with scaled center-of-mass energy z = cme / meta.
pub fn scaled_cross_section_2_4(z: f64) -> (f64, f64, f64) {
    if z <= 4.0 {
        return (0.0, 0.0, 0.0);
    }
    // Scaled masses.
    let masses = vec![1.0, 1.0, 1.0, 1.0];
    let nevents = 500000;
    // Incomping particle momenta.
    let p1 = FourMomentum {
        e: z / 2.0,
        px: 0.0,
        py: 0.0,
        pz: (z * z / 4.0 - 1.0).sqrt(),
    };
    let p2 = FourMomentum {
        e: z / 2.0,
        px: 0.0,
        py: 0.0,
        pz: -(z * z / 4.0 - 1.0).sqrt(),
    };

    let msqrd44 = |fm: &Vec<FourMomentum>| {
        let p3 = fm[0];
        let p4 = fm[1];
        let p5 = fm[2];
        let p6 = fm[3];

        let diag1 = mat_elem_2_to_4_diag_4pt(&p1, &p2, &p3, &p4, &p5, &p6, &(p4 + p5 + p6));
        let diag2 = mat_elem_2_to_4_diag_4pt(&p1, &p2, &p4, &p3, &p5, &p6, &(p3 + p5 + p6));
        let diag3 = mat_elem_2_to_4_diag_4pt(&p1, &p2, &p5, &p4, &p3, &p6, &(p3 + p4 + p6));
        let diag4 = mat_elem_2_to_4_diag_4pt(&p1, &p2, &p6, &p4, &p5, &p3, &(p3 + p4 + p5));
        let diag5 = mat_elem_2_to_4_diag_4pt(&p1, &p3, &p4, &p2, &p5, &p6, &(-p2 + p5 + p6));
        let diag6 = mat_elem_2_to_4_diag_4pt(&p1, &p3, &p5, &p4, &p2, &p6, &(-p2 + p4 + p6));
        let diag7 = mat_elem_2_to_4_diag_4pt(&p1, &p3, &p6, &p4, &p5, &p2, &(-p2 + p4 + p5));
        let diag8 = mat_elem_2_to_4_diag_4pt(&p1, &p4, &p5, &p2, &p3, &p6, &(-p2 + p3 + p6));
        let diag9 = mat_elem_2_to_4_diag_4pt(&p1, &p4, &p6, &p2, &p5, &p3, &(-p2 + p3 + p5));
        let diag10 = mat_elem_2_to_4_diag_4pt(&p1, &p5, &p6, &p4, &p2, &p3, &(-p2 + p3 + p4));

        (diag1 + diag2 + diag3 + diag4 + diag5 + diag6 + diag7 + diag8 + diag9 + diag10).powi(2)
    };
    let msqrd66 = |fm: &Vec<FourMomentum>| {
        let p3 = fm[0];
        let p4 = fm[1];
        let p5 = fm[2];
        let p6 = fm[3];

        mat_elem_2_to_4_diag_6pt(&p1, &p2, &p3, &p4, &p5, &p6).powi(2)
    };
    let msqrd46 = |fm: &Vec<FourMomentum>| {
        let p3 = fm[0];
        let p4 = fm[1];
        let p5 = fm[2];
        let p6 = fm[3];

        let diag1 = mat_elem_2_to_4_diag_4pt(&p1, &p2, &p3, &p4, &p5, &p6, &(p4 + p5 + p6));
        let diag2 = mat_elem_2_to_4_diag_4pt(&p1, &p2, &p4, &p3, &p5, &p6, &(p3 + p5 + p6));
        let diag3 = mat_elem_2_to_4_diag_4pt(&p1, &p2, &p5, &p4, &p3, &p6, &(p3 + p4 + p6));
        let diag4 = mat_elem_2_to_4_diag_4pt(&p1, &p2, &p6, &p4, &p5, &p3, &(p3 + p4 + p5));
        let diag5 = mat_elem_2_to_4_diag_4pt(&p1, &p3, &p4, &p2, &p5, &p6, &(-p2 + p5 + p6));
        let diag6 = mat_elem_2_to_4_diag_4pt(&p1, &p3, &p5, &p4, &p2, &p6, &(-p2 + p4 + p6));
        let diag7 = mat_elem_2_to_4_diag_4pt(&p1, &p3, &p6, &p4, &p5, &p2, &(-p2 + p4 + p5));
        let diag8 = mat_elem_2_to_4_diag_4pt(&p1, &p4, &p5, &p2, &p3, &p6, &(-p2 + p3 + p6));
        let diag9 = mat_elem_2_to_4_diag_4pt(&p1, &p4, &p6, &p2, &p5, &p3, &(-p2 + p3 + p5));
        let diag10 = mat_elem_2_to_4_diag_4pt(&p1, &p5, &p6, &p4, &p2, &p3, &(-p2 + p3 + p4));
        let diag11 = mat_elem_2_to_4_diag_6pt(&p1, &p2, &p3, &p4, &p5, &p6);

        (diag1 + diag2 + diag3 + diag4 + diag5 + diag6 + diag7 + diag8 + diag9 + diag10) * diag11
    };
    let den = 2.0 * z * (z * z - 4.0).sqrt() * 24.0;

    let res44 = Rambo::integrate_phase_space(&masses, z, &msqrd44, nevents);
    let res66 = Rambo::integrate_phase_space(&masses, z, &msqrd66, nevents);
    let res46 = Rambo::integrate_phase_space(&masses, z, &msqrd46, nevents);
    (res44.0 / den, res66.0 / den, res46.0 / den)
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_msqrd() {
        let p1 = FourMomentum {
            e: 1.36347536791516,
            px: 0.48654385793513644,
            py: 0.7439126984517934,
            pz: 0.262552947610851,
        };
        let p2 = FourMomentum {
            e: 1.4464074265120566,
            px: 0.39391413810037945,
            py: 0.5394471294803707,
            pz: 0.803693280903389,
        };
        let p3 = FourMomentum {
            e: 1.432964864561526,
            px: 0.46806713355288165,
            py: 0.011299288058934875,
            pz: 0.9133311489513443,
        };
        let p4 = FourMomentum {
            e: 1.348372605111796,
            px: 0.1998035697146534,
            py: 0.8528997946489925,
            pz: 0.2252757333424067,
        };
        let p5 = FourMomentum {
            e: 1.3799014297797332,
            px: 0.7184764773751884,
            py: 0.5508025975034925,
            pz: 0.2907507625959014,
        };
        let p6 = FourMomentum {
            e: 1.3549774141968551,
            px: 0.7679257948094826,
            py: 0.35155994853026074,
            pz: 0.3502275963416086,
        };

        let diags = vec![
            mat_elem_2_to_4_diag_4pt(&p1, &p2, &p3, &p4, &p5, &p6, &(p4 + p5 + p6)),
            mat_elem_2_to_4_diag_4pt(&p1, &p2, &p4, &p3, &p5, &p6, &(p3 + p5 + p6)),
            mat_elem_2_to_4_diag_4pt(&p1, &p2, &p5, &p4, &p3, &p6, &(p3 + p4 + p6)),
            mat_elem_2_to_4_diag_4pt(&p1, &p2, &p6, &p4, &p5, &p3, &(p3 + p4 + p5)),
            mat_elem_2_to_4_diag_4pt(&p1, &p3, &p4, &p2, &p5, &p6, &(-p2 + p5 + p6)),
            mat_elem_2_to_4_diag_4pt(&p1, &p3, &p5, &p4, &p2, &p6, &(-p2 + p4 + p6)),
            mat_elem_2_to_4_diag_4pt(&p1, &p3, &p6, &p4, &p5, &p2, &(-p2 + p4 + p5)),
            mat_elem_2_to_4_diag_4pt(&p1, &p4, &p5, &p2, &p3, &p6, &(-p2 + p3 + p6)),
            mat_elem_2_to_4_diag_4pt(&p1, &p4, &p6, &p2, &p5, &p3, &(-p2 + p3 + p5)),
            mat_elem_2_to_4_diag_4pt(&p1, &p5, &p6, &p4, &p2, &p3, &(-p2 + p3 + p4)),
            mat_elem_2_to_4_diag_6pt(&p1, &p2, &p3, &p4, &p5, &p6),
        ];

        let mathematica_diags = vec![
            18.05247982398534,
            17.372884288286393,
            18.54833413427864,
            18.80756268121668,
            -24.52857518513328,
            -85.91138791140078,
            -42.74839960071725,
            -131.28985303066247,
            276.4142767501613,
            50.25650487269819,
            -27.175755213202898,
        ];

        for (x, y) in diags.iter().zip(mathematica_diags.iter()) {
            assert!((x - y).abs() < 1e-5);
        }
    }
    #[test]
    fn test_int_msqrd() {
        use ndarray::prelude::*;
        let zs = Array1::<f64>::linspace(4.0, 100.0, 500);
        for z in zs.iter() {
            let int = scaled_cross_section_2_4(*z);
            println!("s, int = {:e}, {:e}", *z, int.0);
        }
    }
}
