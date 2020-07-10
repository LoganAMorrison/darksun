/// Find the root of a function `f` bracketed by (x1, x2) using the
/// Brent-Decker algorithm.
pub fn brent<F>(f: F, x1: f64, x2: f64, tol: f64) -> Result<f64, String>
where
    F: Fn(f64) -> f64,
{
    let itermax = 100;
    let mut a = x1;
    let mut b = x2;
    let mut c = x2;
    let mut d: f64 = 0.0;
    let mut e: f64 = 0.0;
    let mut fa = f(a);
    let mut fb = f(b);
    let mut p: f64;
    let mut q: f64;
    let mut r: f64;
    let mut s: f64;
    let mut tol1: f64;
    let mut xm: f64;

    if fa > 0.0 && fb > 0.0 || fa < 0.0 && fb < 0.0 {
        return Err("Root must be bracketed in brent".to_owned());
    }
    let mut fc = fb;
    for _ in 0..itermax {
        if (fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0) {
            c = a;
            fc = fa;
            e = b - a;
            d = b - a;
        }
        if fc.abs() < fb.abs() {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        tol1 = 2.0 * f64::EPSILON * b.abs() + 0.5 * tol;
        xm = 0.5 * (c - b);
        if xm.abs() <= tol1 || fb == 0.0 {
            return Ok(b);
        }
        if e.abs() >= tol1 && fa.abs() > fb.abs() {
            s = fb / fa;
            if a == c {
                p = 2.0 * xm * s;
                q = 1.0 - s;
            } else {
                q = fa / fc;
                r = fb / fc;
                p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if p > 0.0 {
                q = -q;
            }
            p = p.abs();
            let min1 = 3.0 * xm * q - (tol1 * q).abs();
            let min2 = (e * q).abs();
            if 2.0 * p < min1.min(min2) {
                e = d;
                d = p / q;
            } else {
                d = xm;
                e = d;
            }
        } else {
            d = xm;
            e = d;
        }
        a = b;
        fa = fb;
        if d.abs() > tol1 {
            b += d;
        } else {
            b += tol1.copysign(xm);
        }
        fb = f(b);
    }
    Err("Maximum number of iterations exceeded in brent".to_owned())
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_brent_quad() {
        let f = |x: f64| x - 2.0;
        let sol = brent(f, 0.0, 3.0, 1e-5);
        println!("sol = {:?}", sol);
    }
}
