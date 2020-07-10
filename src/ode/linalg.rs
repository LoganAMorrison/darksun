//! Module for performing decompositions of various matrices.

use ndarray::prelude::*;

/// Perform a triangular decomposition on a real dense matrix `a`. The upper and
/// lower triangular matricies are stored in `a` upon completion and the pivots
/// are stored in `ip`. Returns an error code which will be equal to zero if
/// decomposition is successful and a non-zero number which stage matrix was
/// found to be singular otherwise.
pub(crate) fn dec(n: usize, mut a: ArrayViewMut2<f64>, mut ip: ArrayViewMut1<i32>) -> usize {
    let mut ier = 0;
    ip[n - 1] = 1;
    if n != 1 {
        let nm1 = n - 1;
        for k in 0..nm1 {
            let kp1 = k + 1;
            let mut m = k;
            for i in kp1..n {
                if a[[i, k]].abs() > a[[m, k]].abs() {
                    m = i;
                }
            }
            ip[k] = m as i32;
            let mut t = a[[m, k]];
            if m != k {
                ip[n - 1] = -ip[n - 1];
                a[[m, k]] = a[[k, k]];
                a[[k, k]] = t;
            }
            if t == 0.0 {
                ier = k;
                ip[n - 1] = 0;
                return ier;
            }
            t = t.recip();
            for i in kp1..n {
                a[[i, k]] *= -t;
            }
            for j in kp1..n {
                t = a[[m, j]];
                a[[m, j]] = a[[k, j]];
                a[[k, j]] = t;
                if t != 0.0 {
                    for i in kp1..n {
                        a[[i, j]] += a[[i, k]] * t;
                    }
                }
            }
        }
    }
    if a[[n - 1, n - 1]] == 0.0 {
        ier = n;
        ip[n - 1] = 0;
    }

    return ier;
}

/// Perform a triangular decomposition on a complex matrix with real and
/// imaginary components `ar` and `ai`. The upper and
/// lower triangular matricies are stored in `a` upon completion and the pivots
/// are stored in `ip`. Returns an error code which will be equal to zero if
/// decomposition is successful and a non-zero number which stage matrix was
/// found to be singular otherwise.
pub(crate) fn decc(
    n: usize,
    mut ar: ArrayViewMut2<f64>,
    mut ai: ArrayViewMut2<f64>,
    mut ip: ArrayViewMut1<i32>,
) -> usize {
    let mut ier = 0;
    ip[n - 1] = 1;
    if n != 1 {
        let nm1 = n - 1;
        for k in 0..nm1 {
            let kp1 = k + 1;
            let mut m = k;
            for i in kp1..n {
                if ar[[i, k]].abs() + ai[[i, k]].abs() > ar[[m, k]].abs() + ai[[m, k]].abs() {
                    m = i;
                }
            }
            ip[k] = m as i32;
            let mut tr = ar[[m, k]];
            let mut ti = ai[[m, k]];
            if m != k {
                ip[n - 1] = -ip[n - 1];
                ar[[m, k]] = ar[[k, k]];
                ai[[m, k]] = ai[[k, k]];
                ar[[k, k]] = tr;
                ai[[k, k]] = ti;
            }
            if tr.abs() + ti.abs() == 0.0 {
                ier = k;
                ip[n - 1] = 0;
                return ier;
            }
            let den = tr * tr + ti * ti;
            tr = tr / den;
            ti = -ti / den;
            for i in kp1..n {
                let prodr = ar[[i, k]] * tr - ai[[i, k]] * ti;
                let prodi = ai[[i, k]] * tr + ar[[i, k]] * ti;
                ar[[i, k]] = -prodr;
                ai[[i, k]] = -prodi;
            }
            for j in kp1..n {
                tr = ar[[m, j]];
                ti = ai[[m, j]];
                ar[[m, j]] = ar[[k, j]];
                ai[[m, j]] = ai[[k, j]];
                ar[[k, j]] = tr;
                ai[[k, j]] = ti;
                if tr.abs() + ti.abs() == 0.0 {
                } else if ti == 0.0 {
                    for i in kp1..n {
                        let prodr = ar[[i, k]] * tr;
                        let prodi = ai[[i, k]] * tr;
                        ar[[i, j]] += prodr;
                        ai[[i, j]] += prodi;
                    }
                } else if tr == 0.0 {
                    for i in kp1..n {
                        let prodr = -ai[[i, k]] * ti;
                        let prodi = ar[[i, k]] * ti;
                        ar[[i, j]] += prodr;
                        ai[[i, j]] += prodi;
                    }
                } else {
                    for i in kp1..n {
                        let prodr = ar[[i, k]] * tr - ai[[i, k]] * ti;
                        let prodi = ai[[i, k]] * tr + ar[[i, k]] * ti;
                        ar[[i, j]] += prodr;
                        ai[[i, j]] += prodi;
                    }
                }
            }
        }
    }
    if ar[[n - 1, n - 1]].abs() + ai[[n - 1, n - 1]].abs() == 0.0 {
        ier = n;
        ip[n - 1] = 0;
    }

    return ier;
}

pub(crate) fn sol(n: usize, a: ArrayView2<f64>, mut b: ArrayViewMut1<f64>, ip: ArrayView1<i32>) {
    if n != 1 {
        let nm1 = n - 1;
        for k in 0..nm1 {
            let kp1 = k + 1;
            let m = ip[k];
            let t = b[m as usize];
            b[m as usize] = b[k];
            b[k] = t;
            for i in kp1..n {
                b[i] += a[[i, k]] * t;
            }
        }
        for k in 0..nm1 {
            let km1 = n - k - 2;
            let kb = km1 + 1;
            b[kb] = b[kb] / a[[kb, kb]];
            let t = -b[kb];
            for i in 0..(km1 + 1) {
                b[i] += a[[i, kb]] * t;
            }
        }
    }
    b[0] = b[0] / a[[0, 0]];
}

pub(crate) fn solc(
    n: usize,
    ar: ArrayView2<f64>,
    ai: ArrayView2<f64>,
    mut br: ArrayViewMut1<f64>,
    mut bi: ArrayViewMut1<f64>,
    ip: ArrayView1<i32>,
) {
    if n != 1 {
        let nm1 = n - 1;
        for k in 0..nm1 {
            let kp1 = k + 1;
            let m = ip[k] as usize;
            let tr = br[m];
            let ti = bi[m];
            br[m] = br[k];
            bi[m] = bi[k];
            br[k] = tr;
            bi[k] = ti;
            for i in kp1..n {
                let prodr = ar[[i, k]] * tr - ai[[i, k]] * ti;
                let prodi = ai[[i, k]] * tr + ar[[i, k]] * ti;
                br[i] += prodr;
                bi[i] += prodi;
            }
        }
        for k in 0..nm1 {
            let km1 = n - k - 2;
            let kb = km1 + 1;
            let den = ar[[kb, kb]] * ar[[kb, kb]] + ai[[kb, kb]] * ai[[kb, kb]];
            let mut prodr = br[kb] * ar[[kb, kb]] + bi[kb] * ai[[kb, kb]];
            let mut prodi = bi[kb] * ar[[kb, kb]] - br[kb] * ai[[kb, kb]];
            br[kb] = prodr / den;
            bi[kb] = prodi / den;
            let tr = -br[kb];
            let ti = -bi[kb];
            for i in 0..(km1 + 1) {
                prodr = ar[[i, kb]] * tr - ai[[i, kb]] * ti;
                prodi = ai[[i, kb]] * tr + ar[[i, kb]] * ti;
                br[i] += prodr;
                bi[i] += prodi;
            }
        }
    }
    let den = ar[[0, 0]] * ar[[0, 0]] + ai[[0, 0]] * ai[[0, 0]];
    let prodr = br[0] * ar[[0, 0]] + bi[0] * ai[[0, 0]];
    let prodi = bi[0] * ar[[0, 0]] - br[0] * ai[[0, 0]];
    br[0] = prodr / den;
    bi[0] = prodi / den;
}
